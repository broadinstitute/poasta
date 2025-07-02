# Multi-Piece Gap-Affine Penalties for Partial Order Alignment

## Mathematical Framework

### Current Gap-Affine Model

The current POASTA implementation uses a standard affine gap penalty model:

$$g_{\text{affine}}(k) = \begin{cases} 
0 & \text{if } k = 0 \\
\alpha + \beta \cdot (k-1) & \text{if } k > 0
\end{cases}$$

where:
- $\alpha$ = gap opening penalty
- $\beta$ = gap extension penalty  
- $k$ = gap length

### Generalized Multi-Piece Gap-Affine Model

The generalized n-piece gap-affine model uses multiple breakpoints to create piecewise linear gap penalties:

$$g_{\text{multi-piece}}(k) = \begin{cases}
0 & \text{if } k = 0 \\
\alpha + \sum_{i=1}^{j} \beta_i \cdot \min(k-k_{i-1}, k_i - k_{i-1}) & \text{if } k_{j-1} < k \leq k_j
\end{cases}$$

where:
- $\alpha$ = gap opening penalty
- $\mathbf{k} = [k_0, k_1, k_2, ..., k_{n-1}, k_n]$ = breakpoints with $k_0 = 1$ and $k_n = \infty$
- $\boldsymbol{\beta} = [\beta_1, \beta_2, ..., \beta_n]$ = extension penalties for each piece
- $j$ = piece index such that $k_{j-1} < k \leq k_j$

This simplifies to:
$$g_{\text{multi-piece}}(k) = \alpha + \sum_{i=1}^{j-1} \beta_i \cdot (k_i - k_{i-1}) + \beta_j \cdot (k - k_{j-1})$$

### Special Cases

This framework naturally handles all standard models as special cases:

**1-piece (Standard Affine)**: $\mathbf{k} = [1, \infty]$, $\boldsymbol{\beta} = [\beta]$
**2-piece**: $\mathbf{k} = [1, k_0, \infty]$, $\boldsymbol{\beta} = [\beta_1, \beta_2]$  
**3-piece**: $\mathbf{k} = [1, k_0, k_1, \infty]$, $\boldsymbol{\beta} = [\beta_1, \beta_2, \beta_3]$

Typically, $\beta_1 > \beta_2 > \beta_3 > ...$ to progressively reduce penalty for longer gaps.

### Biological Motivation

Research has shown that gap length distributions exhibit complex patterns:
- **Short gaps (1-3)**: High penalty per base (structural constraints)
- **Medium gaps (4-20)**: Moderate penalty (common evolutionary events)  
- **Long gaps (>20)**: Low penalty (large structural rearrangements)

Multi-piece models can capture these biological realities:
- **2-piece**: Distinguish short vs long gaps
- **3-piece**: Add ultra-long gap category for large structural variants
- **n-piece**: Model complex empirical distributions with arbitrary precision

## Algorithmic Design for Graph Alignment

### State Space Modification

Current POASTA uses a 3-state model:
- $M(i,v)$ = Match state at query position $i$, graph node $v$
- $I(i,v)$ = Insertion state (gap in graph)  
- $D(i,v)$ = Deletion state (gap in query)

For multi-piece gap-affine, we need to track both gap lengths and piece indices:

#### Generalized State Space

$$S = \{M(i,v), I(i,v,j,l), D(i,v,j,l)\}$$

where:
- $M(i,v)$ = Match state
- $I(i,v,j,l)$ = Insertion state, gap in piece $j$ of length $l$
- $D(i,v,j,l)$ = Deletion state, gap in piece $j$ of length $l$
- $j \in [1, n]$ = piece index 
- $l$ = gap length within current piece

This unified representation handles any number of pieces without separate state types.

### Dynamic Programming Recurrence Relations

#### Match State
$$M(i,v) = \min_{u \in \text{pred}(v)} \left\{ M(i-1,u) + s(q_i, v), \min_{j,l} \{I(i-1,u,j,l), D(i-1,u,j,l)\} + s(q_i, v) \right\}$$

where $s(q_i, v)$ is the substitution score.

#### Insertion States

**Starting new gap:**
$$I(i,v,1,1) = M(i-1,v) + \alpha + \beta_1$$

**Extending within same piece:**
$$I(i,v,j,l+1) = I(i-1,v,j,l) + \beta_j \quad \text{if } l+1 \leq k_j$$

**Transitioning to next piece:**
$$I(i,v,j+1,k_{j-1}+1) = I(i-1,v,j,k_{j-1}) + \beta_{j+1} \quad \text{if } j < n$$

#### Deletion States

**Starting new gap:**
$$D(i,v,1,1) = \min_{u \in \text{pred}(v)} M(i,u) + \alpha + \beta_1$$

**Extending within same piece:**
$$D(i,v,j,l+1) = \min_{u \in \text{pred}(v)} D(i,u,j,l) + \beta_j \quad \text{if } l+1 \leq k_j$$

**Transitioning to next piece:**
$$D(i,v,j+1,k_{j-1}+1) = \min_{u \in \text{pred}(v)} D(i,u,j,k_{j-1}) + \beta_{j+1} \quad \text{if } j < n$$

## Implementation Strategy

### 1. Data Structure Changes

#### Generalized VisitedCell Structure
```rust
struct VisitedCellMultiPieceAffine {
    visited_m: Score,
    visited_i: Vec<Vec<Score>>,  // visited_i[piece][length]
    visited_d: Vec<Vec<Score>>,  // visited_d[piece][length]
}
```

#### Flexible Gap Scoring Structure
```rust
pub struct GapMultiPieceAffine {
    cost_mismatch: u8,
    cost_gap_open: u8,
    breakpoints: Vec<usize>,     // [k₁, k₂, ..., kₙ₋₁] 
    extension_costs: Vec<u8>,    // [β₁, β₂, ..., βₙ]
}

impl GapMultiPieceAffine {
    // Factory methods for common cases
    pub fn standard_affine(open: u8, extend: u8) -> Self {
        Self {
            cost_mismatch: 4,
            cost_gap_open: open,
            breakpoints: vec![],
            extension_costs: vec![extend],
        }
    }
    
    pub fn two_piece(open: u8, short_extend: u8, long_extend: u8, breakpoint: usize) -> Self {
        Self {
            cost_mismatch: 4,
            cost_gap_open: open,
            breakpoints: vec![breakpoint],
            extension_costs: vec![short_extend, long_extend],
        }
    }
    
    pub fn three_piece(open: u8, costs: [u8; 3], breakpoints: [usize; 2]) -> Self {
        Self {
            cost_mismatch: 4,
            cost_gap_open: open,
            breakpoints: breakpoints.to_vec(),
            extension_costs: costs.to_vec(),
        }
    }
}
```

### 2. Modified A* Algorithm

#### Generalized State Representation
```rust
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AlignmentStateMultiPiece {
    pub query_offset: usize,
    pub graph_node: POANodeIndex,
    pub state_type: AlignmentStateType,
    pub piece_index: usize,  // Which piece we're in
    pub gap_length: usize,   // Length within current piece
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AlignmentStateType {
    Match,
    Insertion,
    Deletion,
}
```

#### Generalized Gap Cost Calculation
```rust
impl GapMultiPieceAffine {
    fn gap_cost(&self, gap_length: usize) -> usize {
        if gap_length == 0 {
            return 0;
        }
        
        let mut cost = self.cost_gap_open as usize;
        let mut remaining = gap_length;
        let mut prev_breakpoint = 1;
        
        for (i, &breakpoint) in self.breakpoints.iter().enumerate() {
            let piece_length = std::cmp::min(remaining, breakpoint - prev_breakpoint + 1);
            cost += self.extension_costs[i] as usize * piece_length;
            remaining -= piece_length;
            prev_breakpoint = breakpoint + 1;
            
            if remaining == 0 {
                break;
            }
        }
        
        // Handle remaining length in final piece
        if remaining > 0 {
            let final_piece_idx = self.extension_costs.len() - 1;
            cost += self.extension_costs[final_piece_idx] as usize * remaining;
        }
        
        cost
    }
    
    fn find_piece(&self, gap_length: usize) -> usize {
        for (i, &breakpoint) in self.breakpoints.iter().enumerate() {
            if gap_length <= breakpoint {
                return i;
            }
        }
        self.extension_costs.len() - 1  // Final piece
    }
}
```

### 3. Generalized Heuristic Function

The heuristic function works with any number of pieces:

```rust
fn min_gap_cost_multi_piece(&self, remaining_query: usize) -> usize {
    if remaining_query == 0 {
        return 0;
    }
    
    // Use the gap cost calculation directly - it's already optimal
    self.gap_model.gap_cost(remaining_query)
}
```

### 4. Memory Optimization

To manage the increased memory requirements:

1. **Lazy State Creation**: Only create gap states when needed
2. **Gap Length Pruning**: Limit maximum tracked gap length
3. **Blocked Storage**: Extend current blocked storage to handle multiple gap states
4. **State Compression**: Use bit packing for small gap lengths

### 5. Parameter Selection Examples

**Standard Affine (1-piece):**
```rust
GapMultiPieceAffine::standard_affine(6, 2)
```

**Two-piece (recommended):**
```rust
GapMultiPieceAffine::two_piece(6, 2, 1, 3)
// α=6, β₁=2, β₂=1, k₀=3
```

**Three-piece (for ultra-long gaps):**
```rust
GapMultiPieceAffine::three_piece(6, [2, 1, 0], [3, 20])
// Short gaps (1-3): penalty 2
// Medium gaps (4-20): penalty 1  
// Long gaps (>20): penalty 0 (free extension)
```

## Complexity Analysis

### Time Complexity
- Current: $O(|V| \cdot |E| \cdot n)$ 
- Multi-piece: $O(|V| \cdot |E| \cdot n \cdot k_{\max} \cdot p)$

where $k_{\max}$ is the maximum gap length tracked and $p$ is the number of pieces.

### Space Complexity  
- Current: $O(|V| \cdot n)$
- Multi-piece: $O(|V| \cdot n \cdot k_{\max} \cdot p)$

The complexity scales linearly with the number of pieces, making it practical for small $p$ (2-4 pieces).

## Implementation Plan

### Phase 1: Core Data Structures
1. Implement `GapMultiPieceAffine` scoring model with factory methods
2. Create generalized `VisitedCellMultiPieceAffine` structure
3. Extend `AlignmentState` to track piece indices and gap lengths

### Phase 2: Algorithm Modification
1. Update A* recurrence relations for multi-piece transitions
2. Implement generalized heuristic functions
3. Add piece transition logic

### Phase 3: Optimization
1. Add memory optimization for sparse piece states
2. Implement dynamic piece selection based on gap length
3. Performance testing across different piece configurations

### Phase 4: Integration
1. Update configuration system with multi-piece support
2. Add parameter validation and sensible defaults
3. Update CLI interface and documentation

## Key Advantages of Generalized Design

1. **Unified Framework**: All gap penalty models (linear, affine, multi-piece) handled by same code
2. **Easy Configuration**: Factory methods make common cases simple
3. **Extensible**: Adding new pieces requires no algorithm changes
4. **Backward Compatible**: Standard affine is a special case
5. **Research Friendly**: Easy to experiment with different piece configurations

This generalized design provides maximum flexibility while maintaining code simplicity and mathematical rigor.