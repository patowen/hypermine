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

/*
Proof that simple nodes strings are well behaved:

Notation:
* A or `A` === string or substring (uppercase) (Backticks used to resolve potential ambiguity)
* a or `a` === generator (lowercase) (Backticks used to resolve potential ambiguity)
* @A or @(A) === simplified version of A based on "append" rules. "@" has precedence over concatenation
* A+a === A after appending a using simplification rules. Note: @(Aa) = @A+a. Concatenation has precedence over "+"
* A=B === A and B are the same string (same generators in the same order)
* A!=B === not (A=B)
* a<b === a has priority over b. (Justification: view priority as a rank) Note, if a<b, then a~b, as priority is not defined for non-adjacent characters
* a<=b === a<b or a=b
* a~b === a and b commute (including equality). Note that this is NOT an equivalence relation.
* a~A === a commutes with everything in A

Definition of "append":
    AaB + a = AB iff `a` strongly slots into AaB before aB
    AC + a = AaC iff `a` strongly slots into AC before C, and C!=aB for any B

Definition of "strongly slots": `a` strongly slots into AB before B iff
    * `a` weakly slots into AB before B, and
    * If `a` weakly slots into CDB before DB, then D is empty. (That is, it doesn't weakly slot any earlier)

Definition of "weakly slots": `a` weakly slots into AB before B iff
    * a~B, and
    * B != bC with b<a (Note: If B is nonempty, this is equivalent to B = bC with a<=b)

Lemma 1: The "append" rules are well-defined (as in, there's exactly one option)
Proof of lemma:
    For every string, there is exactly one location in any string A that some generator `a` can strongly slot into.
    Proof:
        `a` weakly slots into any string A before "" (as it, it slots in at the end) because the definition of "weakly slots" holds trivially
        Since "strongly slots" just means "is the leftmost spot where it weakly slots", we've got the minimum of a finite nonzero set of indices, which is unique.
    Therefore, if we split A into A=BC, where `a` strongly slots into A before C, then C is defined uniquely.
    C either starts with `a`, or it doesn't. (An empty C is seen as not starting with `a`). Both of these cases reach one branch of the definition of "append".
    Therefore, "append" is well defined, so it's not an abuse of notation to say A+a=B for some A, a, and B.

Definition: A is simple iff A does not contain the forbidden substring bAa where a<=b and a~b and a~A

Lemma 2: The "append" rules turn any simple string into some simple string
Proof of lemma:
    Assume A is simple. WTS: A+a is simple for all `a`
    Case 1: A = BaC with a~C
        Then, C does not contain `a` because that would form a forbidden substring in A: aDa with a~D
        `a` weakly slots into A before aC because a~C and a<=a
        `a` doesn't weakly slot in any earlier
        Proof:
            Let DE = B be defined such that `a` weakly slots into A (=DEaC) before EaC.
            Suppose, by way of contradiction, that E is nonempty. That is, let E=eF
                Then, a~e and a~F and a<=e
                This gives us the forbidden substring of `A`: eFa, which violates our initial constraint that `A` is simple.
                By contradiciton, E must be empty.
        `a` strongly slots into A before aC
        Therefore, A+a = BC
        Removing a character from a simple string keeps it simple, since no new opportunities arise for forbidden substrings (proof is relatively trivial and will be skipped)
        Therefore, in Case 1, A+a is simple.
    Case 2: The first case doesn't hold.
        Let BC = A be defined such that `a` strongly slots into A before C.
        C does not contain `a` because a~C, so this would give us case 1, which we constrained to not hod for case 2.
        Therefore, A+a = BaC by definition of "append".
        To check that BaC is simple, we just need to check that there's no b in B or c in C that can create a forbidden substring with `a`.
        Proof for b in B: Suppose, by way of contradiction, that B ends in bD where bDa is a forbidden substring.
            Then, a<=b, a~b, and a~D. Note that we also have a~C from before.
            Therefore, since A ends in bDC, and a~bDC, and a<=b this means that `a` weakly slots into A before bDC.
            This contradicts the fact that `a` strongly slots into A before C. This completes the proof by contradiction.
        Proof for c in C: Suppose, by way of contradiction, that C=DcE where aDc is a forbidden substring.
            Then, c<=a and c~a and c~D. Case 1 doesn't hold, so c!=a, so c<a.
            D is nonempty because otherwise, this would contradict the fact that `a` weakly slots into A before C (=DcE = cE).
            So, if we define dF = D, then `a` weakly slots into A before dFcE. This means that a<=d by definition of "weakly slots"
            Since A is simple, dFc cannot be a forbidden substring. We have c~F and c~d because c~D, so this means that d<c (because (not c<=d) and d~c).
            a<=d and d<c and a~c, so by the transitive property of priority, a<c. This contradicts the fact that c<a. This completes the proof by contradiction.
        Therefore, in Case 2, A+a is simple.
    Therefore, A+a is simple for all `a`, completing the proof of lemma 2.

TODO: What still needs to be proven:
    1. The "append" rules do not change the node being pointed to. That is, A+a is equivalent to Aa in a group-theoretic sense. This is quite easy to show.
    2. No two distinct simple strings are equivalent in a group-theoretic sense. By the regional division geometric intuition, this can likely be accepted without proof, but a proof would be good to be more certain.
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

            if new_segment < segment {
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

            if new_segment < segment {
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
    #[ignore] // This test is not a real test, as it is just used to discover a formula for how many nodes are "x" steps away with a given node string prefix.
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
        let mut graph = Graph::<NodeString>::default();
        *graph.get_mut(NodeId::ROOT) = Some(NodeString::default());
        let mut ns = NodeString::new();
        test_no_overlap_recursive(&mut graph, NodeId::ROOT, &mut ns, 7);
    }

    fn test_no_overlap_recursive(
        g: &mut Graph<NodeString>,
        node: NodeId,
        ns: &mut NodeString,
        len: u32,
    ) {
        if len == 0 {
            return;
        }

        for side in Side::iter() {
            if ns.has_child(side) {
                if let Some(old_neighbor) = g.neighbor(node, side) {
                    assert!(
                        g.get(old_neighbor).is_none(), /* Assert we haven't seen this node yet. It's fine if it has already been created as a side-effect of ensure_neighbor */
                        "Reached a graph node twice: {:?} vs {:?}",
                        ns,
                        g.get(old_neighbor)
                    );
                }

                ns.append(side);

                let neighbor = g.ensure_neighbor(node, side);
                *g.get_mut(neighbor) = Some(ns.clone());
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
