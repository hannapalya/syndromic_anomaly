# Basic supervised LSTM + validation/test split + metrics + bank holidays
# - No z-scores / AE / GMM / k-of-m / deseasonalizing
# - Global scaler fit ONLY on TRAIN sims' first TRAIN_DAYS
# - Train on FULL length of train sims (labels: outbreak anywhere in window)
# - Validation/Test on last 49 weeks of holdout sims (labels: endpoint)
# - Metrics match your R functions; also bank-holiday-masked versions

import os, numpy as np, pandas as pd, tensorflow as tf
from sklearn.metrics import confusion_matrix

# ===== CONFIG =====
DATA_DIR      = "signal_datasets"
BANK_HOLS_CSV = "Bankholidays.csv"  # must contain a 'date' column (YYYY-MM-DD) or first col is dates
SIGNALS       = list(range(1,17))
DAYS_PER_YEAR = 364
SEQ, STRIDE   = 14, 1
TRAIN_YEARS   = 6
TRAIN_DAYS    = TRAIN_YEARS * DAYS_PER_YEAR
VALID_DAYS    = 49 * 7
EPOCHS        = 30
BATCH_SIZE    = 128
RNG_STATE     = 42
OVERSAMPLE_POS= 3
START_DATE    = "2000-01-01"  # used if your CSVs lack explicit dates

# ===== HELPERS =====
def load_data(sig):
    X = pd.read_csv(os.path.join(DATA_DIR, f"simulated_totals_sig{sig}.csv"))
    Y = (pd.read_csv(os.path.join(DATA_DIR, f"simulated_outbreaks_sig{sig}.csv")) > 0).astype(int)
    # optional date column support
    date_col = None
    for c in ["date", "Date", "ds", "timestamp"]:
        if c in X.columns: date_col = c; break
    dates = pd.to_datetime(X[date_col]).dt.normalize().to_numpy() if date_col else None
    if date_col:  # drop date col from X/Y if present
        X = X.drop(columns=[date_col])
        if date_col in Y.columns: Y = Y.drop(columns=[date_col])
    return X, Y, dates

def stratified_sim_split(sims, rng, train_frac=0.5, hold_val_frac=0.5):
    pos_idx = [i for i,s in enumerate(sims) if s["y"].sum() > 0]
    neg_idx = [i for i,s in enumerate(sims) if s["y"].sum() == 0]
    rng.shuffle(pos_idx); rng.shuffle(neg_idx)
    n_train_pos = int(len(pos_idx) * train_frac)
    n_train_neg = int(len(neg_idx) * train_frac)
    train_ids = pos_idx[:n_train_pos] + neg_idx[:n_train_neg]
    hold_ids  = pos_idx[n_train_pos:] + neg_idx[n_train_neg:]
    rng.shuffle(hold_ids)
    n_val = int(len(hold_ids) * (1 - hold_val_frac))
    val_ids  = hold_ids[:n_val]
    test_ids = hold_ids[n_val:]
    return ([sims[i] for i in train_ids],
            [sims[i] for i in val_ids],
            [sims[i] for i in test_ids])

def make_seq_labels(x, y, start, end, seq_len=SEQ, stride=STRIDE, label_mode="within_window"):
    xs, ys = [], []
    for i in range(start, end - seq_len + 1, stride):
        xs.append(x[i:i+seq_len])
        if label_mode == "within_window":
            ys.append(int(y[i:i+seq_len].max() > 0))
        else:  # endpoint
            ys.append(int(y[i+seq_len-1] > 0))
    if not xs:
        return np.empty((0, seq_len, 1), np.float32), np.empty((0,), np.int32)
    return np.asarray(xs, np.float32)[..., None], np.asarray(ys, np.int32)

def build_lstm(seq_len=SEQ):
    inp = tf.keras.Input(shape=(seq_len,1))
    x = tf.keras.layers.LSTM(64)(inp)
    x = tf.keras.layers.Dense(32, activation='relu')(x)
    out = tf.keras.layers.Dense(1, activation='sigmoid')(x)
    m = tf.keras.Model(inp, out)
    m.compile(optimizer="adam", loss="binary_crossentropy",
              metrics=[tf.keras.metrics.AUC(name="auc")])
    return m

def sens_spec(y_true, y_pred):
    TN, FP, FN, TP = confusion_matrix(y_true, y_pred, labels=[0,1]).ravel()
    sens = TP/(TP+FN) if (TP+FN)>0 else 0.0
    spec = TN/(TN+FP) if (TN+FP)>0 else 0.0
    return sens, spec

def tune_threshold(y_val, p_val):  # Youden's J
    best_t, best_j = 0.5, -1
    for t in np.linspace(0.05, 0.95, 19):
        yhat = (p_val >= t).astype(int)
        s, sp = sens_spec(y_val, yhat)
        j = s + sp - 1
        if j > best_j:
            best_j, best_t = j, float(t)
    return best_t

def oversample(X, Y, k=OVERSAMPLE_POS):
    if k <= 1: return X, Y
    pos = (Y == 1)
    if pos.sum() == 0: return X, Y
    Xp, Yp = X[pos], Y[pos]
    X_aug = np.concatenate([X] + [Xp]*(k-1), axis=0)
    Y_aug = np.concatenate([Y] + [Yp]*(k-1), axis=0)
    idx = np.random.permutation(len(Y_aug))
    return X_aug[idx], Y_aug[idx]

# --- R-equivalent metrics (A: T×J alarms, O: day labels aligned to A via tail) ---
def _align_tail(O, T): return O[-T:, :] if O.shape[0] >= T else O

def compute_specificity(A, O):
    Oa = _align_tail(O, A.shape[0])
    TN = np.logical_and(A == 0, Oa == 0).sum()
    FP = np.logical_and(A == 1, Oa == 0).sum()
    return (TN / (TN + FP)) if (TN + FP) > 0 else np.nan

def compute_sensitivity(A, O):
    Oa = _align_tail(O, A.shape[0])
    TP = np.logical_and(A == 1, Oa > 0).sum()
    FN = np.logical_and(A == 0, Oa > 0).sum()
    return (TP / (TP + FN)) if (TP + FN) > 0 else np.nan

def compute_pod(A, O):
    Oa = _align_tail(O, A.shape[0])
    return np.mean((np.logical_and(A == 1, Oa > 0)).sum(axis=0) > 0)

def compute_timeliness(A, O):
    Oa = _align_tail(O, A.shape[0])
    T, J = A.shape
    score = 0.0
    for j in range(J):
        y = Oa[:, j]; a = A[:, j]
        idx_out = np.where(y > 0)[0]
        if len(idx_out) == 0: score += 1.0; continue
        idx_hit = np.where((a == 1) & (y > 0))[0]
        if len(idx_hit) == 0: score += 1.0; continue
        r1, r2 = int(idx_out[0]), int(idx_out[-1])
        obs = int(idx_hit[0])
        score += (obs - r1) / (r2 - r1 + 1)
    return score / J

# --- Bank holiday support ---
def load_bank_holidays(path=BANK_HOLS_CSV):
    df = pd.read_csv(path)
    col = "date" if "date" in df.columns else df.columns[0]
    return set(pd.to_datetime(df[col]).dt.normalize())

def build_holiday_mask_for_tail(total_len, valid_days, seq, bank_set, dates_series=None, start_date=START_DATE):
    if dates_series is not None:
        all_dates = pd.to_datetime(dates_series).dt.normalize().to_numpy()
    else:
        all_dates = pd.date_range(start=start_date, periods=total_len, freq="D").normalize().to_numpy()
    tail_dates = all_dates[-valid_days:]
    end_dates = tail_dates[seq-1:]  # window ends
    return np.array([d in bank_set for d in end_dates], dtype=bool)

def compute_specificity_masked(A, O, M):
    Oa = _align_tail(O, A.shape[0]); keep = ~M
    TN = np.logical_and.reduce([A == 0, Oa == 0, keep]).sum()
    FP = np.logical_and.reduce([A == 1, Oa == 0, keep]).sum()
    return (TN / (TN + FP)) if (TN + FP) > 0 else np.nan

def compute_sensitivity_masked(A, O, M):
    Oa = _align_tail(O, A.shape[0]); keep = ~M
    TP = np.logical_and.reduce([A == 1, Oa > 0, keep]).sum()
    FN = np.logical_and.reduce([A == 0, Oa > 0, keep]).sum()
    return (TP / (TP + FN)) if (TP + FN) > 0 else np.nan

def compute_pod_masked(A, O, M):
    Oa = _align_tail(O, A.shape[0]); any_hit = []
    for j in range(A.shape[1]):
        keep = ~M[:, j]
        any_hit.append((np.logical_and.reduce([A[:, j] == 1, Oa[:, j] > 0, keep])).any())
    return np.mean(any_hit)

def compute_timeliness_masked(A, O, M):
    Oa = _align_tail(O, A.shape[0]); J = A.shape[1]; score = 0.0
    for j in range(J):
        keep = ~M[:, j]; y = Oa[:, j][keep]; a = A[:, j][keep]
        idx_out = np.where(y > 0)[0]
        if len(idx_out) == 0: score += 1.0; continue
        idx_hit = np.where((a == 1) & (y > 0))[0]
        if len(idx_hit) == 0: score += 1.0; continue
        r1, r2 = int(idx_out[0]), int(idx_out[-1]); obs = int(idx_hit[0])
        score += (obs - r1) / (r2 - r1 + 1)
    return score / J

# ===== MAIN =====
np.random.seed(RNG_STATE); tf.random.set_seed(RNG_STATE)
rng = np.random.RandomState(RNG_STATE)
bank_set = load_bank_holidays(BANK_HOLS_CSV)

summary_all, summary_bh = {}, {}

for S in SIGNALS:
    print(f"\n--- Signal {S} ---")
    Xsig, Ysig, dates_all = load_data(S)

    sims = []
    for sim_idx, col in enumerate(Xsig.columns):
        x = Xsig[col].to_numpy(float)
        y = Ysig[col].to_numpy(int)
        if len(x) >= TRAIN_DAYS + VALID_DAYS:
            sims.append(dict(x=x, y=y, sim=f"sig{S}_sim{sim_idx}",
                             dates=dates_all if dates_all is not None else None))
    if not sims:
        print("No complete sims; skip."); continue

    train_sims, val_sims, test_sims = stratified_sim_split(sims, rng, train_frac=0.5, hold_val_frac=0.5)

    # Global scaler from TRAIN sims' first TRAIN_DAYS
    all_train_hist = np.concatenate([s["x"][:TRAIN_DAYS] for s in train_sims])
    g_mn, g_mx = float(all_train_hist.min()), float(all_train_hist.max())
    g_den = (g_mx - g_mn) if g_mx > g_mn else 1.0
    def scale(x): return ((x - g_mn) / g_den).astype(np.float32)

    # TRAIN (full length, within-window labels), with oversampling
    XtrL, YtrL = [], []
    for d in train_sims:
        x = scale(d["x"]); y = d["y"]
        Xt, Yt = make_seq_labels(x, y, 0, len(x), SEQ, STRIDE, label_mode="within_window")
        if len(Xt): XtrL.append(Xt); YtrL.append(Yt)
    if not XtrL: print("No training sequences."); continue
    Xtr = np.concatenate(XtrL); Ytr = np.concatenate(YtrL)
    Xtr, Ytr = oversample(Xtr, Ytr, OVERSAMPLE_POS)

    # VAL: last 49w, endpoint labels
    XvL, YvL = [], []
    for d in val_sims:
        x = scale(d["x"]); y = d["y"]; s = len(x) - VALID_DAYS
        Xv, Yv = make_seq_labels(x, y, s, len(x), SEQ, STRIDE, label_mode="endpoint")
        if len(Xv): XvL.append(Xv); YvL.append(Yv)
    Xval = np.concatenate(XvL) if XvL else np.empty((0,SEQ,1), np.float32)
    Yval = np.concatenate(YvL) if YvL else np.empty((0,), np.int32)

    # TEST per-sim (keep separate for metrics)
    per_sim_preds, per_sim_labels, per_sim_masks = [], [], []
    Xte_concat, Yte_concat = [], []
    for d in test_sims:
        x = scale(d["x"]); y = d["y"]; s = len(x) - VALID_DAYS
        Xt, Yt = make_seq_labels(x, y, s, len(x), SEQ, STRIDE, label_mode="endpoint")
        per_sim_labels.append(Yt.astype(int))
        Xte_concat.append(Xt); Yte_concat.append(Yt)
        # build BH mask per sim for window ends
        M = build_holiday_mask_for_tail(total_len=len(x), valid_days=VALID_DAYS, seq=SEQ,
                                        bank_set=bank_set, dates_series=d["dates"])
        per_sim_masks.append(M.astype(bool))
    Xte = np.concatenate(Xte_concat) if Xte_concat else np.empty((0,SEQ,1), np.float32)
    Yte = np.concatenate(Yte_concat) if Yte_concat else np.empty((0,), np.int32)

    # Fit
    tf.keras.backend.clear_session()
    model = build_lstm(SEQ)
    es = tf.keras.callbacks.EarlyStopping(monitor="val_loss" if len(Xval) else "loss",
                                          patience=5, restore_best_weights=True)
    if len(Xval):
        model.fit(Xtr, Ytr, validation_data=(Xval, Yval), epochs=EPOCHS,
                  batch_size=BATCH_SIZE, verbose=0, callbacks=[es])
        Pval = model.predict(Xval, verbose=0).squeeze()
        t_star = tune_threshold(Yval, Pval)
    else:
        model.fit(Xtr, Ytr, epochs=EPOCHS, batch_size=BATCH_SIZE, verbose=0, callbacks=[es])
        t_star = 0.5
    print(f"Chosen threshold t* = {t_star:.2f}")

    # Predict test (both concatenated and per-sim)
    if len(Xte):
        Pte = model.predict(Xte, verbose=0).squeeze()
        yhat_concat = (Pte >= t_star).astype(int)

        # split back to per-sim preds
        ofs = 0; per_sim_preds.clear()
        for Yt in per_sim_labels:
            n = len(Yt); per_sim_preds.append(yhat_concat[ofs:ofs+n]); ofs += n

        # Build A (T×J), O (T×J), M (T×J)
        A = np.stack(per_sim_preds, axis=1)
        O = np.stack(per_sim_labels, axis=1)
        M = np.stack(per_sim_masks, axis=1)

        # Scalar metrics (all days)
        sens = compute_sensitivity(A, O)
        spec = compute_specificity(A, O)
        pod  = compute_pod(A, O)
        tim  = compute_timeliness(A, O)
        print(f"ALL DAYS → Sens={sens:.3f}, Spec={spec:.3f}, POD={pod:.3f}, Tim={tim:.3f}")

        # Masked metrics (exclude bank holidays)
        sens_bh = compute_sensitivity_masked(A, O, M)
        spec_bh = compute_specificity_masked(A, O, M)
        pod_bh  = compute_pod_masked(A, O, M)
        tim_bh  = compute_timeliness_masked(A, O, M)
        print(f"NON-BANK-HOLIDAYS → Sens={sens_bh:.3f}, Spec={spec_bh:.3f}, POD={pod_bh:.3f}, Tim={tim_bh:.3f}")

        summary_all[S] = dict(sensitivity=sens, specificity=spec, pod=pod, timeliness=tim, thr=t_star)
        summary_bh[S]  = dict(sensitivity=sens_bh, specificity=spec_bh, pod=pod_bh, timeliness=tim_bh, thr=t_star)
    else:
        print("No test sequences.")

# ===== SUMMARY =====
if summary_all:
    df = pd.DataFrame.from_dict(summary_all, orient="index")
    print("\n=== SUMMARY (ALL DAYS) ==="); print(df); print("\nMeans:\n", df.mean(numeric_only=True))
if summary_bh:
    dfb = pd.DataFrame.from_dict(summary_bh, orient="index")
    print("\n=== SUMMARY (NON-BANK-HOLIDAYS) ==="); print(dfb); print("\nMeans:\n", dfb.mean(numeric_only=True))

# ===== SUMMARY =====
if summary_all:
    df = pd.DataFrame.from_dict(summary_all, orient="index")
    print("\n=== SUMMARY (ALL DAYS) ==="); print(df); print("\nMeans:\n", df.mean(numeric_only=True))
    
    # Save results to CSV
    df.to_csv("LSTM_results_all_days.csv")
    print(f"\nResults saved to: LSTM_results_all_days.csv")
    
if summary_bh:
    dfb = pd.DataFrame.from_dict(summary_bh, orient="index")
    print("\n=== SUMMARY (NON-BANK-HOLIDAYS) ==="); print(dfb); print("\nMeans:\n", dfb.mean(numeric_only=True))
    
    # Save results to CSV
    dfb.to_csv("LSTM_results_non_bank_holidays.csv")
    print(f"Results saved to: LSTM_results_non_bank_holidays.csv")

# Create a combined summary with both metrics
if summary_all and summary_bh:
    combined_summary = {}
    for signal in summary_all.keys():
        combined_summary[f"Signal_{signal}"] = {
            "Sensitivity_All": summary_all[signal]["sensitivity"],
            "Specificity_All": summary_all[signal]["specificity"],
            "POD_All": summary_all[signal]["pod"],
            "Timeliness_All": summary_all[signal]["timeliness"],
            "Threshold_All": summary_all[signal]["thr"],
            "Sensitivity_NonBH": summary_bh[signal]["sensitivity"],
            "Specificity_NonBH": summary_bh[signal]["specificity"],
            "POD_NonBH": summary_bh[signal]["pod"],
            "Timeliness_NonBH": summary_bh[signal]["timeliness"],
            "Threshold_NonBH": summary_bh[signal]["thr"]
        }
    
    df_combined = pd.DataFrame.from_dict(combined_summary, orient="index")
    df_combined.to_csv("LSTM_results_first.csv")
    print(f"Combined results saved to: LSTM_results_first.csv")
    
    # Print overall performance summary
    print(f"\n=== OVERALL PERFORMANCE SUMMARY ===")
    print(f"Average Sensitivity (All Days): {df['sensitivity'].mean():.3f}")
    print(f"Average Specificity (All Days): {df['specificity'].mean():.3f}")
    print(f"Average POD (All Days): {df['pod'].mean():.3f}")
    print(f"Average Timeliness (All Days): {df['timeliness'].mean():.3f}")
    print(f"Average Sensitivity (Non-BH): {dfb['sensitivity'].mean():.3f}")
    print(f"Average Specificity (Non-BH): {dfb['specificity'].mean():.3f}")
    print(f"Average POD (Non-BH): {dfb['pod'].mean():.3f}")
    print(f"Average Timeliness (Non-BH): {dfb['timeliness'].mean():.3f}")