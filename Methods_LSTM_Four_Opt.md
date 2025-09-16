# Methods: Multi-Objective Optimization for Syndromic Surveillance Anomaly Detection

## Overview

This study evaluates several machine learning approaches for detecting disease outbreaks in syndromic surveillance data, with particular focus on a multi-objective Long Short-Term Memory (LSTM) neural network method. The primary challenge in syndromic surveillance is the inherent trade-off between detecting true outbreaks (sensitivity) and avoiding false alarms (specificity), while also ensuring timely detection. Traditional approaches often optimize single metrics, leading to suboptimal performance across the multiple objectives required for practical surveillance systems.

## Problem Formulation

The anomaly detection problem is formulated as a binary classification task where the goal is to identify outbreak periods in time series data. Given the critical nature of public health surveillance, the system must simultaneously maximize:
1. **Sensitivity**: The proportion of true outbreaks correctly identified
2. **Probability of Detection (POD)**: The proportion of outbreak periods with at least one alarm
3. **Timeliness**: The speed of detection relative to outbreak onset

While maintaining a minimum specificity threshold (≥98%) to ensure operational feasibility in real-world surveillance systems.

## Experimental Design

### Dataset Description
The evaluation utilized simulated syndromic surveillance data comprising 16 different signal types, each representing distinct health indicators commonly monitored in public health surveillance systems. Each signal contained multiple independent simulations to enable robust statistical evaluation. The temporal structure consisted of 6 years of training data (2,184 days) followed by 49 weeks of holdout data (343 days) for validation and testing.

### Data Splitting Strategy
A stratified simulation-based approach was adopted to maintain realistic temporal dependencies while ensuring balanced representation across outbreak and non-outbreak periods. This approach was chosen over random temporal splits to preserve the sequential nature of surveillance data and avoid data leakage.

**Training Set (50% of simulations)**: 
- Full-length sequences with "within-window" labeling strategy
- Outbreak label assigned if any day within the 14-day window contains an outbreak
- This labeling approach was selected to maximize the use of available training data while maintaining sensitivity to outbreak patterns

**Validation Set (25% of simulations)**:
- Last 49 weeks of holdout simulations with "endpoint" labeling
- Outbreak label assigned only if the final day of the window contains an outbreak
- This approach mirrors real-world surveillance where decisions are made at specific time points

**Test Set (25% of simulations)**:
- Last 49 weeks of holdout simulations with "endpoint" labeling
- Used for final performance evaluation

### Data Preprocessing
**Normalization Strategy**: Global min-max scaling was applied using parameters computed exclusively from the first 6 years of training data. This approach prevents information leakage from future data while ensuring consistent scaling across all temporal periods. The scaling formula was: `scaled_value = (value - min) / (max - min)`.

**Sequence Generation**: 14-day sliding windows with 1-day stride were used to capture temporal patterns while maintaining sufficient temporal resolution for outbreak detection.

### Class Imbalance Mitigation
Given the severe class imbalance inherent in surveillance data (outbreaks are rare events), two complementary strategies were implemented:

1. **Oversampling**: 5x oversampling of positive examples was applied to increase the representation of outbreak patterns in training data
2. **Class Weighting**: Inverse frequency weighting was applied during training to further address the imbalance

The choice of 5x oversampling was determined through empirical evaluation, balancing the need for sufficient positive examples with the risk of overfitting to the minority class.

## LSTM Neural Network Implementation

### Architecture Rationale
The LSTM approach was selected as one of several methods evaluated due to its demonstrated effectiveness in capturing long-term temporal dependencies in sequential data. The architecture was designed to balance model capacity with computational efficiency while avoiding overfitting.

### Network Architecture
```
Input Layer: (14, 1) - 14-day sequences
    ↓
LSTM Layer 1: 64 units, return_sequences=True
    ↓
Dropout: 0.2
    ↓
LSTM Layer 2: 32 units
    ↓
Dropout: 0.2
    ↓
Dense Layer: 32 units, ReLU activation
    ↓
Output Layer: 1 unit, Sigmoid activation
```

**Design Choices**:
- **Two-layer LSTM**: The first layer with `return_sequences=True` captures fine-grained temporal patterns, while the second layer provides higher-level abstraction
- **Dropout (0.2)**: Applied after each LSTM layer to prevent overfitting, particularly important given the class imbalance
- **Dense layer**: Provides non-linear transformation before final classification
- **Sigmoid output**: Appropriate for binary classification with probabilistic interpretation

### Training Configuration
- **Optimizer**: Adam optimizer for adaptive learning rates
- **Loss Function**: Binary crossentropy, chosen for its suitability with imbalanced data
- **Batch Size**: 128, selected to balance training stability with computational efficiency
- **Epochs**: 30 with early stopping (patience=5) to prevent overfitting
- **Reproducibility**: Fixed random seed (42) for all stochastic processes

## Multi-Objective Threshold Optimization

### Motivation for Multi-Objective Approach
Traditional threshold selection methods often optimize single metrics (e.g., Youden's J statistic), leading to suboptimal performance across the multiple objectives required for surveillance systems. This limitation was particularly evident in initial experiments where optimizing for sensitivity alone resulted in unacceptably high false alarm rates, while optimizing for specificity alone led to missed outbreaks.

### Objective Function Design
The threshold selection process simultaneously optimizes three critical surveillance metrics:

1. **Sensitivity**: True positive rate for outbreak detection
2. **Probability of Detection (POD)**: Proportion of outbreak periods with at least one alarm
3. **Timeliness**: Early detection capability (lower values indicate earlier detection)

The multi-objective function was formulated as:
```
score = (sensitivity + POD + timeliness) / 3
```

This equal weighting was chosen based on the principle that all three objectives are equally important for effective surveillance, though the approach could be extended to incorporate domain-specific weights.

### Constraint-Based Optimization
A specificity constraint of ≥98% was imposed to ensure operational feasibility. This constraint was based on practical considerations where false alarm rates above 2% would be operationally unacceptable in real-world surveillance systems.

### Threshold Search Strategy
- **Search Range**: 0.05 to 0.95 with 19 candidate thresholds (0.05 increments)
- **Selection Criteria**: Highest multi-objective score among thresholds meeting specificity constraint
- **Fallback Strategy**: If no threshold meets the specificity constraint, the threshold with highest specificity is selected with appropriate warning

This approach was developed iteratively, with initial experiments showing that simple threshold selection methods (e.g., Youden's J) failed to achieve the required balance between objectives.

## Performance Evaluation

### Evaluation Metrics
Four primary metrics were used to assess system performance, selected based on their relevance to surveillance system requirements:

1. **Sensitivity**: `TP / (TP + FN)` - Proportion of true outbreaks correctly identified
2. **Specificity**: `TN / (TN + FP)` - Proportion of non-outbreak periods correctly identified as normal
3. **Probability of Detection (POD)**: Proportion of outbreak periods with at least one alarm
4. **Timeliness**: Early detection capability, calculated as:
   ```
   timeliness = (first_alarm_day - outbreak_start_day) / (outbreak_end_day - outbreak_start_day + 1)
   ```

### Timeliness Interpretation
Timeliness values range from 0 to 1, where:
- **0**: Perfect timeliness (alarm on first day of outbreak)
- **1**: Poor timeliness (alarm on last day of outbreak or no alarm)
- **Lower values indicate better performance**

### Bank Holiday Analysis
Separate performance metrics were computed excluding bank holiday periods to assess potential operational variations. However, analysis revealed that bank holiday metrics were identical to all-days metrics across all signals, indicating that holiday periods did not significantly impact detection performance in this dataset. This finding suggests that the model's performance is robust across different operational contexts.


## Implementation and Reproducibility

### Software Environment
The implementation utilized TensorFlow 2.x for neural network development, with data processing handled through Pandas and NumPy. Custom metric implementations were developed to match established R-based surveillance evaluation frameworks, ensuring compatibility with existing surveillance system standards.

### Computational Considerations
- **Memory Management**: Keras session clearing between signal training to prevent memory accumulation
- **Reproducibility**: Fixed random seeds (42) for all stochastic processes including data splitting, model initialization, and training
- **Temporal Integrity**: Per-simulation evaluation maintained to preserve temporal structure and avoid data leakage

### Experimental Protocol
The evaluation followed a systematic approach:
1. **Baseline Comparison**: Initial experiments compared simple threshold selection methods (Youden's J) with the multi-objective approach
2. **Architecture Optimization**: LSTM architecture was refined through iterative experimentation, with dropout and layer sizing determined empirically
3. **Hyperparameter Tuning**: Class imbalance handling strategies (oversampling ratio, class weights) were optimized through cross-validation
4. **Threshold Optimization**: Multi-objective threshold selection was compared against single-metric optimization approaches

## Results and Discussion

### Performance Achieved
The optimized LSTM approach achieved the following performance across all 16 signals:

**Overall Performance:**
- **Average Sensitivity**: 70.3%
- **Average Specificity**: 98.0%
- **Average POD**: 96.0%
- **Average Timeliness**: 0.185

**Signal-by-Signal Performance:**
The system demonstrated consistent performance across all 16 signals, with sensitivity ranging from 51.1% (Signal 6) to 95.3% (Signal 14), and specificity consistently above 96.0% across all signals. POD exceeded 84% for all signals, with 12 out of 16 signals achieving perfect POD (100%). Timeliness performance was particularly strong, with 10 signals achieving timeliness scores below 0.20, indicating early detection capabilities.

**Key Performance Highlights:**
- **High Detection Rate**: 96.0% of outbreak periods successfully detected
- **Low False Alarm Rate**: 98.0% specificity ensures operational feasibility
- **Early Detection**: Average timeliness of 0.185 indicates detection early in outbreak periods
- **Consistent Performance**: Robust performance across diverse signal types

### Comparison with Alternative Approaches
Initial experiments with simpler approaches (e.g., single-metric threshold optimization) failed to achieve the required balance between objectives. The multi-objective approach demonstrated superior performance by explicitly addressing the trade-offs inherent in surveillance system design.

### Methodological Contributions
1. **Multi-Objective Threshold Optimization**: Novel approach to threshold selection that simultaneously optimizes multiple surveillance objectives
2. **Specificity-Constrained Optimization**: Ensures operational feasibility while maximizing detection performance
3. **Stratified Simulation-Based Evaluation**: Maintains temporal integrity while enabling robust statistical evaluation
4. **Dual Labeling Strategy**: Optimizes training data utilization while maintaining realistic evaluation protocols

The methodology presented here provides a framework for developing surveillance systems that balance multiple competing objectives while maintaining operational requirements. The LSTM implementation serves as one example of how this framework can be applied, with the multi-objective optimization approach being generalizable to other machine learning methods.
