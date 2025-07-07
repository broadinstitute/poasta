use tracing_subscriber::filter::FilterFn;

pub fn align_state_filter() -> FilterFn {
    FilterFn::new(|metadata| {
        let is_poasta = metadata.target().starts_with("poasta::aligner");
        let is_visit = metadata.target().ends_with("set_visited");
        let is_extend = metadata.target().ends_with("extend");

        if is_poasta && (is_visit || is_extend) {
            return false;
        }
        
        true
    })
}

pub fn only_align_states() -> FilterFn {
    FilterFn::new(|metadata| {
        let is_poasta = metadata.target().starts_with("poasta::aligner");
        if !is_poasta {
            return false;
        }
        
        let is_visit = metadata.target().ends_with("set_visited");
        let is_queue = metadata.target().ends_with("queue_item");
        let is_extend = metadata.target().ends_with("extend");
        let is_new_seq = metadata.name() == "align_seq";
        
        is_visit || is_queue || is_new_seq || is_extend
    })
}