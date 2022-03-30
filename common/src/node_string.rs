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
* A+a === A after appending `a` using simplification rules. Note: @(Aa) = @A+a. Concatenation has precedence over "+"
* A+B === A after appending all generators in B one at a time using simplification rules. Note: @(AB) = @A+B. Concatenation has precedence over "+"
* A=B === A and B are the same string (same generators in the same order)
* A!=B === not (A=B)
* a<b === a has priority over b. (Justification: view priority as a rank) Note, if a<b, then a$b, as priority is not defined for non-commuting characters
* a<=b === a<b or a=b
* a$b === a and b commute (including equality). Note that this is NOT an equivalence relation, but it is symmetric.
* a$A === a commutes with everything in A. Elements within A do NOT have to mutually commute.
* A~B === A and B correspond to the same group element in the relevant right-angled Coxeter group.
    Note that this is an equivalence relation, so rules related to equivalence relations will be taken as given.
    We also take group theory rules as given, such as associativity, identity (``), and inverse, so some steps will be skipped based on that.

Definition of "append":
    AaB + a = AB iff `a` strongly slots into AaB before aB
    AC + a = AaC iff `a` strongly slots into AC before C, and C!=aB for any B

Definition of "strongly slots": `a` strongly slots into AB before B iff
    * `a` weakly slots into AB before B, and
    * If `a` weakly slots into CDB before DB, then D is empty. (That is, it doesn't weakly slot any earlier)

Definition of "weakly slots": `a` weakly slots into AB before B iff
    * a$B, and
    * B != bC with b<a (Note: If B is nonempty, this is equivalent to B = bC with a<=b)

Definition of the "@" symbol (Note: the notation for this signifies simplifying a node string)
    * @(``) = `` (that is, the base case of this recursive definition is that simplifying an empty string gives the empty string)
    * @(Aa) = @A+a (this part of the definition is stating that elements are appended one at a time)

Definition of "~" (equivalence in Coxeter group)
    * AB ~ AaaB (Reflecting about the same element twice does nothing)
    * If a$b, then AB ~ AababB (Going all the way around a vertex does nothing)
    * For any given A, you can find any B~A with a finite number of applications of the above rules (This proof is not formal, so this is just to allow an induction proof to work.)

Definition: A is simple iff A does not contain the forbidden substring bAa where a<=b and a$b and a$A

Lemma: The "append" rules are well-defined (as in, there's exactly one option)
Proof of lemma:
    For every string, there is exactly one location in any string A that some generator `a` can strongly slot into.
    Proof:
        `a` weakly slots into any string A before "" (as it, it slots in at the end) because the definition of "weakly slots" holds trivially
        Since "strongly slots" just means "is the leftmost spot where it weakly slots", we've got the minimum of a finite nonzero set of indices, which is unique.
    Therefore, if we split A into A=BC, where `a` strongly slots into A before C, then C is defined uniquely.
    C either starts with `a`, or it doesn't. (An empty C is seen as not starting with `a`). Both of these cases reach one branch of the definition of "append".
    Therefore, "append" is well defined, so it's not an abuse of notation to say A+a=B for some A, a, and B.

Simplification Lemma: The "append" rules turn any simple string into some simple string (This also shows that for any A, @A is simple, and this won't be explicitly proven.)
Proof of lemma:
    Assume A is simple. WTS: A+a is simple for all `a`
    Case 1: A = BaC with a$C
        Then, C does not contain `a` because that would form a forbidden substring in A: aDa with a$D
        `a` weakly slots into A before aC because a$C and a<=a
        `a` doesn't weakly slot in any earlier
        Proof:
            Let DE = B be defined such that `a` weakly slots into A (=DEaC) before EaC.
            Suppose, by way of contradiction, that E is nonempty. That is, let E=eF
                Then, a$e and a$F and a<=e
                This gives us the forbidden substring of `A`: eFa, which violates our initial constraint that `A` is simple.
                By contradiciton, E must be empty.
        `a` strongly slots into A before aC
        Therefore, A+a = BC
        Removing a character from a simple string keeps it simple, since no new opportunities arise for forbidden substrings (proof is relatively trivial and will be skipped)
        Therefore, in Case 1, A+a is simple.
    Case 2: The first case doesn't hold.
        Let BC = A be defined such that `a` strongly slots into A before C.
        C does not contain `a` because a$C, so this would give us case 1, which we constrained to not hod for case 2.
        Therefore, A+a = BaC by definition of "append".
        To check that BaC is simple, we just need to check that there's no b in B or c in C that can create a forbidden substring with `a`.
        Proof for b in B: Suppose, by way of contradiction, that B ends in bD where bDa is a forbidden substring.
            Then, a<=b, a$b, and a$D. Note that we also have a$C from before.
            Therefore, since A ends in bDC, and a$bDC, and a<=b this means that `a` weakly slots into A before bDC.
            This contradicts the fact that `a` strongly slots into A before C. This completes the proof by contradiction.
        Proof for c in C: Suppose, by way of contradiction, that C=DcE where aDc is a forbidden substring.
            Then, c<=a and c$a and c$D. Case 1 doesn't hold, so c!=a, so c<a.
            D is nonempty because otherwise, this would contradict the fact that `a` weakly slots into A before C (=DcE = cE).
            So, if we define dF = D, then `a` weakly slots into A before dFcE. This means that a<=d by definition of "weakly slots"
            Since A is simple, dFc cannot be a forbidden substring. We have c$F and c$d because c$D, so this means that d<c (because (not c<=d) and d$c).
            a<=d and d<c and a$c, so by the transitive property of priority, a<c. This contradicts the fact that c<a. This completes the proof by contradiction.
        Therefore, in Case 2, A+a is simple.
    Therefore, A+a is simple for all `a`, completing the proof of the lemma.

Idempotence of Simplification Lemma: If A is simple, @A = A. A direct corollary of this lemma is that @@A = @A.
Proof of lemma:
     We just need to show that if Aa is simple, then A+a = Aa. The lemma will follow by induction with a trivial base case.
     To do this, we show that `a` can only weakly slot into A at the very end. (It follows immediately that `a` can also only strongly slot into A at the very end, forcing the append rules to yield A+a = Aa)
        Suppose, by way of contradiction, that A=BbC, and `a` weakly slots into BbC before bC. By definition, a$bC and a<=b. Since a$bC, that means that a$b and a$C.
        However, Aa=BbCa, and bCa is a forbidden substring because a<=b and a$b and a$C. This contradicts the fact that Aa is simple, proving the lemma.

Commutivity Lemma: If a$A, then aA ~ Aa
Proof of lemma:
    To show that a$b implies ab ~ ba, use the following equivalence chain: ab ~ a(abab)b ~ (aa)ba(aa) ~ ba
    This can be used as the inductive step for the full lemma. The base case (empty A) is trivial, so we'll skip it.
    Inductive step: Assume true for A and prove for Ab: To be explicit, the assumption is that a$A implies aA ~ Aa
        To prove for Ab, start with a$Ab. In other words, a$A and a$b. Then, we have aAb ~ Aab (by assumption) ~ Aba (by the initial thing we showed). This completes the inductive step, and the proof.

Correctness of Simplification Lemma: The "append" rules do not change the node being pointed to. That is, A+a ~ Aa
Proof of lemma:
    Case 1: Suppose we have AaB + a = AB with `a` strongly slotting into AaB before aB. In this case, we want to show that AB ~ AaBa
        Then, `a` also weakly slots into AaB before aB, which means that a$B. By the commutivity lemma, then, we have aB ~ Ba.
        This gives AaBa ~ AaaB ~ AB, completing the proof of this case.
    Case 2: Suppose we have AC + a = AaC with `a` strongly slotting into AC before C. In this case, we want to show that AaC ~ ACa
        With the exact same reasoning as before, we have aC ~ Ca. This gives us AaC ~ ACa immediately, completing the full proof.

Strong Slot Stability Lemma: If `a` strongly slots into AB before B, and `a` weakly slots into AC before C, then `a` strongly slots into AC before C.
Proof of lemma:
    Let A = DbE. We just need to show that `a` doesn't weakly slot any earlier into AC, that is, `a` doesn't weakly slot into DbEC before bEC.
    Suppose by way of contradiction that `a` weakly slots into DbEC before bEC. Then a<=b and a$bEC.
    However, since `a` strongly slots into AB before B, we have a$B which implies a$bEB.
    Since a<=b and a$bEB, that means `a` weakly slots into DbEB (=AB) before bEB, contradicting the fact that `a` strongly slots into AB before B. This completes the proof.

Simple Slot Lemma: If AaB is a simple string, and a$B, then `a` strongly slots into AB before B.
Proof of lemma:
    First, we want to show that `a` weakly slots into AB before B.
        Suppose by way of contradiction that `a` does not weakly slot into AB before B.
        Since we know that a$B, this means that the other part of the "weakly slot" definition must fail.
        In other words, B = bC with b<a. However, this would mean that AaB = AabC, creating the forbidden string `ab` with b<=a and b$a. This shows that `a` must weakly slot into AB before B.
    Next, we show that `a` cannot weakly slot any earlier.
        Suppose by way of contradiction that we can: That is, if we let CbD = A, then `a` weakly slots into CbDB before bDB.
        This would mean that a$bDB and a<=b. But then, AaB = CbDaB would have the forbidden string bDa, making AaB not simple. This shows that `a` must strongly slot into AB before B, completing the proof.

Uniqueness of Simplification Lemma: If A~B, then @A = @B.
Proof outline:
    The proof takes advantage of the fact for all A, B, the following holds: @(AB) = @(@(A)B). This can be shown via induction with the idempotence of simplification lemma.
    To prove this, we need to show two things.
        1. @(Aaa) = @A for all A, a
        2. @(Aabab) = @A for all A and for all a$b
    For an example of why this is sufficient, (1) shows that @(AaaB) = @(@(Aaa)B) = @(@(A)B) = @(AB), so having `aa` inserted at the end does not sacrifice generality.
    These two things can actually be reduced to two simpler things by applying earlier lemmas:
        1. A+aa = A for all A, a where A is simple
        2. A+abab = A for all A, a, b where a$b and a!=b and A is simple (restricting to a!=b is allowed because otherwise, we're covered by A+aa = A)
    Based on the recursive definition of "~", an inductive proof would be relatively trivial and is omitted for brevity.
Proof of lemma:
    Part 1: Proof that A+aa = A for all A, a where A is simple.
        Case 1: Suppose A = BaC, where BaC + a = BC (first part of the append definition). We want to show that BC + a = BaC = A.
            We know by the "append" definition that `a` strongly slots into BaC before aC. It directly follows that a$C. Since BaC is simple, that shows by the simple slot lemma that `a` strongly slots into BC before C.
            C does not start with A, because that would create a forbidden substring `aa` in A. Hence, BC + a = BaC = A, finishing case 1.
        Case 2: Suppose A = BC, and BC + a = BaC (second part of the append definition). We want to show that BaC + a = BC = A
            We know by the "append" definition that `a` strongly slots into BC before C. It directly follows that a$C.
            As a$a and a<=a, we immediately have that `a` weakly slots into BaC before aC. By the strong slot stability lemma, `a` strongly slots into BaC before aC, proving that BaC + a = BC = A.
        These two cases complete part 1 of the proof.
    Part 2: Proof that A+abab = A for all A, a, b where a$b and a!=b and A is simple (This part is a lot more informal because the approach changed last-minute to avoid an 8-case proof)
        To prove this, we introduce the concept of a generator's slot. TODO: See if the rest of the proof can be refactored with this concept whether this can make this proof more formal without being more verbose.
        The main idea is that any simple string A has a unique slot for each generator (e.g. `a` and `b`). This slot can be filled or empty.
            The slot for a generator `a` in A is filled if `a` strongly slots into A before a substring that starts with `a`. Otherwise, it is empty.
        Part 1 of the proof of the uniqueness of simplification lemma effectively demonstrated that repeatedly appending a generator toggles its slot between full and empty but doesn't move it.
        To prove Part 2, we want to show that appending `b`, which toggles `b`'s slot, doesn't move `a`'s slot.
        There is potential ambiguity if `a`'s slot and `b`'s slot are both empty: Which is before the other? For the purpose of this proof, the higher-priority generator is first.
        Case 1: `a`'s slot is before `b`'s slot.
            Since a$b, and the strong slot stability lemma holds, unless the two slots are adjacent, there is nothing changing that can cause `a`'s slot to move.
            If the two slots are adjacent, then in all scenarios, we're forced to have a<b due to the fact that `a`'s slot is before `b`'s slot.
                By the definition of "weakly slots" and simple string, there cannot be a generator `c` after `b`'s slot with c<b
                If `b` toggles from full to empty, then `a`'s slot does not move because `a` still weakly slots to the same location, as any possible generator `c` after `a`'s slot would have b<=c, so a<=c (because a<=b)
                If `b` toggles from empty to full, then `a`'s slot does not move because `a` weakly slots before `b`.
        Case 2: `a`'s slot is after `b`'s slot.
            If `b` toggles from full to empty, there is no reason for `a`'s slot to move because a$b, which means that `b` wasn't blocking `a` from moving earlier.
            If `b` toggles from empty to full, this could only potentially cause `a`'s slot to move if a<b, in which case `a`'s slot would have the potential to move right before `b`. Consider this case:
                The two slots cannot be adjacent, since we would be forced to have b<a, which would not cause `a`'s slot to move.
                In this case `a` commutes with `b` and everything after it. One of these generators is right after `b` but before `a`. Call it `c`.
                We have b<=c (by definition of weakly slots for `b`), and we have a<b, so by transitivity, a<c.
                However, since `a` commutes with `c` and everything after it, and a<c, `a`'s original slot would have been incorrect, since there would have been an earlier valid slot.
                This contradiction means that `a`'s slot can't have moved after all.
        These two cases cover everything.
        If we append with A+abab, all that will do is toggle the `a` slot and the `b` slot twice without moving either of them. This means that A+abab = A

Final theorem: @A = @B if and only if A~B
This follows directly from the uniqueness of simplification lemma and the correctness of simplification lemma
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
                ns.append(side);

                if let Some(old_neighbor) = g.neighbor(node, side) {
                    assert!(
                        g.get(old_neighbor).is_none(), /* Assert we haven't seen this node yet. It's fine if it has already been created as a side-effect of ensure_neighbor */
                        "Reached a graph node twice: {:?} vs {:?}",
                        ns.path,
                        g.get(old_neighbor).as_ref().unwrap().path
                    );
                }

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
