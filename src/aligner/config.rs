use super::cost_models::AlignmentCostModel;

pub trait AlignerConfig {
    type CostModel: AlignmentCostModel;
}
