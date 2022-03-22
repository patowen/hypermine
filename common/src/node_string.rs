use serde::{Deserialize, Serialize};

use crate::dodeca::Side;

/* Notes:
Number of child nodes at a given number of steps follows these recurrence relations:
Contribution of face: [n+1] = 9*([n] - [n-1]) + [n-2] with start (1), 1, 6, 46 (or 8*[n] - [n-1] - 1)
Contribution of edge: [n+1] = 9*([n] - [n-1]) + [n-2] with start 0, 1, 8 (or 8*[n] - [n-1])
Contribution of vert: [n+1] = 9*([n] - [n-1]) + [n-2] with start 0, 0, 1 (or 8*[n] - [n-1] + 1)

Matrix:
   0   1   0
   0   0   1
   1  -9   9

Eigs:
L^3 - 9L^2 + 9L - 1 = 0
(L - 1)(L^2 - 8L + 1) = 0
L = 1, 4 +- sqrt(15) (Let p = 4 + sqrt(15))

Explicit formula: [n] = a + bp^n + cp^(-n)

For face:
a = 1/6, b = 5/12 + 1/12*sqrt(15), c = 5/12 - 1/12*sqrt(15)

For edge:
a = 0, b = 1/30*sqrt(15), c = -1/30*sqrt(15)

For vertex:
a = -1/6, b = 1/12 - 1/60*sqrt(15), c = 1/12 + 1/60*sqrt(15)
*/

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct NodeString {
    path: Vec<Side>,
}

impl NodeString {
    pub fn new() -> Self {
        Self { path: vec![] }
    }

    pub fn append(&mut self, new_segment: Side) {
        let mut insertion_point = self.path.len();

        for (index, &segment) in self.path.iter().enumerate().rev() {
            if segment == new_segment {
                self.path.remove(index);
                return;
            }

            if !segment.adjacent_to(new_segment) {
                break;
            }

            if new_segment > segment {
                insertion_point = index;
            }
        }

        self.path.insert(insertion_point, new_segment);
    }

    pub fn backtrack(&mut self) {
        self.path.pop();
    }

    pub fn has_child(&self, new_segment: Side) -> bool {
        for &segment in self.path.iter().rev() {
            if segment == new_segment {
                return false;
            }

            if !segment.adjacent_to(new_segment) {
                return true;
            }

            if new_segment > segment {
                return false;
            }
        }
        true
    }

    pub fn len(&self) -> usize {
        self.path.len()
    }

    pub fn is_empty(&self) -> bool {
        self.path.is_empty()
    }
}

impl Default for NodeString {
    fn default() -> Self {
        NodeString::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        dodeca::Side,
        graph::{Graph, NodeId},
    };

    use super::NodeString;

    #[test]
    fn vertex_sides() {
        let mut ns = NodeString::new();
        // ns.append(Side::B);
        // ns.append(Side::A);
        for i in 0..10 {
            std::println!("{}", num_leaves(&mut ns, i));
        }
    }

    fn num_leaves(ns: &mut NodeString, len: u32) -> u32 {
        if len == 0 {
            return 1;
        }

        let mut count = 0;
        for side in Side::iter() {
            if ns.has_child(side) {
                ns.append(side);
                count += num_leaves(ns, len - 1);
                ns.backtrack();
            }
        }
        count
    }

    #[test]
    fn test_no_overlap() {
        let mut graph = Graph::<()>::default();
        let mut ns = NodeString::new();
        test_no_overlap_recursive(&mut graph, NodeId::ROOT, &mut ns, 9);
    }

    fn test_no_overlap_recursive(g: &mut Graph<()>, node: NodeId, ns: &mut NodeString, len: u32) {
        if len == 0 {
            return;
        }

        for side in Side::iter() {
            if ns.has_child(side) {
                assert!(
                    g.neighbor(node, side).is_none(),
                    "Reached a graph node twice: {:?}",
                    ns
                );

                let neighbor = g.ensure_neighbor(node, side);

                ns.append(side);
                test_no_overlap_recursive(g, neighbor, ns, len - 1);
                ns.backtrack();
            }
        }
    }

    #[test]
    fn test_append_simplifies() {
        let ns = NodeString::new();
        test_append_simplifies_recursive(&ns, 7)
    }

    fn test_append_simplifies_recursive(ns: &NodeString, len: u32) {
        if len == 0 {
            return;
        }

        for side in Side::iter() {
            let mut ns_clone = ns.clone();
            ns_clone.append(side);
            assert_simple(&ns_clone);
            test_append_simplifies_recursive(&ns_clone, len - 1);
        }
    }

    fn assert_simple(ns: &NodeString) {
        for i in 0..ns.len() {
            let test_ns = NodeString {
                path: ns.path[0..i].to_vec(),
            };
            assert!(
                test_ns.has_child(ns.path[i]),
                "Node string not simple: {:?}",
                ns
            );
        }
    }
}
