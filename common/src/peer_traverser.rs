use std::sync::LazyLock;

use crate::{
    array_vec::ArrayVec,
    dodeca::Side,
    graph::{Graph, NodeId},
};

pub struct PeerNode {
    node_id: NodeId,
    parent_path: ArrayVec<Side, 2>,
    child_path: ArrayVec<Side, 2>,
}

impl PeerNode {
    #[inline]
    pub fn node(&self) -> NodeId {
        self.node_id
    }

    #[inline]
    pub fn path_from_peer(&self) -> impl ExactSizeIterator<Item = Side> + use<> {
        self.parent_path.into_iter().rev()
    }

    #[inline]
    pub fn path_from_base(&self) -> impl ExactSizeIterator<Item = Side> + use<> {
        self.child_path.into_iter()
    }
}

pub struct PeerTraverser {
    parent_path: ArrayVec<Side, 2>,
    parent_path_nodes: ArrayVec<NodeId, 3>,
    child_path_index: usize,
}

impl PeerTraverser {
    pub fn new(base_node: NodeId) -> Self {
        let mut parent_path_nodes = ArrayVec::new_with_default(NodeId::ROOT);
        parent_path_nodes.push(base_node);
        PeerTraverser {
            parent_path: ArrayVec::new_with_default(Side::A),
            parent_path_nodes,
            child_path_index: 0,
        }
    }

    /// Assuming `parent_path` obeys shortlex rules up to right before the last
    /// element for a given depth,
    fn parent_path_end_is_shortlex(&self, depth: usize) -> bool {
        if depth <= 1 {
            // One-element node strings are always shortlex.
            return true;
        };
        let last = self.parent_path[depth - 1];
        let second_last = self.parent_path[depth - 2];
        if last == second_last {
            // Backtracking is not valid (short part of shortlex)
            return false;
        }
        if last.adjacent_to(second_last) && (last as usize) < (second_last as usize) {
            // Unnecessarily having a higher side index first is not valid (lex part of shortlex)
            return false;
        }
        true
    }

    /// Assuming `parent_path` and `parent_path_nodes` is already valid apart from possibly the last node,
    /// iterates to the next valid path for the given depth.
    /// Returns `false` if this is not possible.
    /// If `allow_unchanged_path` is true, the `parent_path` will not be incremented if it is already valid, but
    /// the `parent_path_nodes` still will.
    #[must_use]
    fn increment_parent_path_for_depth(
        &mut self,
        graph: &Graph,
        depth: usize,
        mut allow_unchanged_path: bool,
    ) -> bool {
        if depth == 0 {
            // Empty paths are always valid, but they cannot be incremented.
            return allow_unchanged_path;
        }
        loop {
            if allow_unchanged_path && self.parent_path_end_is_shortlex(depth) {
                if let Some(node) = graph.neighbor(
                    self.parent_path_nodes[depth - 1],
                    self.parent_path[depth - 1],
                ) {
                    if graph.length(self.parent_path_nodes[depth - 1]) == graph.length(node) + 1 {
                        self.parent_path_nodes[depth] = node;
                        return true;
                    }
                }
            }

            let mut current_side = self.parent_path[depth - 1];
            current_side = Side::VALUES[(current_side as usize + 1) % Side::VALUES.len()]; // Cycle the current side
            if current_side == Side::A {
                // We looped, so make sure to increment an earlier part of the path.
                if !self.increment_parent_path_for_depth(graph, depth - 1, false) {
                    return false;
                }
            }
            self.parent_path[depth - 1] = current_side;

            allow_unchanged_path = true; // The path has changed, so it won't necessarily need to be changed again.
        }
    }

    #[must_use]
    fn increment_parent_path(&mut self, graph: &Graph) -> bool {
        if self.parent_path.len() == 0 {
            self.parent_path.push(Side::A);
            self.parent_path_nodes.push(NodeId::ROOT);
            return self.increment_parent_path_for_depth(graph, 1, true);
        } else if self.parent_path.len() == 1 {
            if self.increment_parent_path_for_depth(graph, 1, false) {
                return true;
            }
            self.parent_path.fill(Side::A);
            self.parent_path.push(Side::A);
            self.parent_path_nodes.push(NodeId::ROOT);
            return self.increment_parent_path_for_depth(graph, 1, true)
                && self.increment_parent_path_for_depth(graph, 2, true);
        } else if self.parent_path.len() == 2 {
            return self.increment_parent_path_for_depth(graph, 2, false);
        }
        false
    }

    #[must_use]
    fn increment_child_path(
        &mut self,
        graph: &mut impl GraphRef,
        mut allow_unchanged_path: bool,
    ) -> Option<PeerNode> {
        if self.parent_path.len() == 1 {
            let child_paths = &DEPTH1_CHILD_PATHS[self.parent_path[0] as usize];
            loop {
                if allow_unchanged_path {
                    if self.child_path_index >= child_paths.len() {
                        return None;
                    }
                    let child_side = child_paths[self.child_path_index];
                    let mut current_node = self.parent_path_nodes[1];
                    current_node = graph.neighbor(current_node, child_side);
                    if graph.length(current_node) == graph.length(self.parent_path_nodes[0]) {
                        let mut result_child_path = ArrayVec::new_with_default(Side::A);
                        result_child_path.push(child_side);
                        return Some(PeerNode {
                            node_id: current_node,
                            parent_path: self.parent_path,
                            child_path: result_child_path,
                        });
                    }
                }
                self.child_path_index += 1;
                allow_unchanged_path = true;
            }
        } else if self.parent_path.len() == 2 {
            let child_paths =
                &DEPTH2_CHILD_PATHS[self.parent_path[0] as usize][self.parent_path[1] as usize];
            loop {
                if allow_unchanged_path {
                    if self.child_path_index >= child_paths.len() {
                        return None;
                    }
                    let child_path = &child_paths[self.child_path_index];
                    let mut current_node = self.parent_path_nodes[2];
                    for &side in child_path {
                        current_node = graph.neighbor(current_node, side); // TODO
                    }
                    if graph.length(current_node) == graph.length(self.parent_path_nodes[0]) {
                        let mut result_child_path = ArrayVec::new_with_default(Side::A);
                        result_child_path.push(child_path[0]);
                        result_child_path.push(child_path[1]);
                        return Some(PeerNode {
                            node_id: current_node,
                            parent_path: self.parent_path,
                            child_path: result_child_path,
                        });
                    }
                }
                self.child_path_index += 1;
                allow_unchanged_path = true;
            }
        }
        None
    }

    fn next_impl(&mut self, mut graph: impl GraphRef) -> Option<PeerNode> {
        let mut allow_unchanged_path = false;
        loop {
            if let Some(node) = self.increment_child_path(&mut graph, allow_unchanged_path) {
                return Some(node);
            }
            if !self.increment_parent_path(graph.as_ref()) {
                return None;
            }
            allow_unchanged_path = true;
            self.child_path_index = 0;
        }
    }

    pub fn next(&mut self, graph: &Graph) -> Option<PeerNode> {
        self.next_impl(ImmutableGraphRef { graph })
    }

    pub fn ensure_next(&mut self, graph: &mut Graph) -> Option<PeerNode> {
        self.next_impl(MutableGraphRef { graph })
    }
}

static DEPTH1_CHILD_PATHS: LazyLock<[ArrayVec<Side, 5>; Side::VALUES.len()]> =
    LazyLock::new(|| {
        Side::VALUES.map(|parent_side| {
            let mut path_list: ArrayVec<Side, 5> = ArrayVec::new_with_default(Side::A);
            for child_side in Side::iter() {
                if !child_side.adjacent_to(parent_side) {
                    continue;
                }
                path_list.push(child_side);
            }
            path_list
        })
    });

static DEPTH2_CHILD_PATHS: LazyLock<
    [[ArrayVec<[Side; 2], 2>; Side::VALUES.len()]; Side::VALUES.len()],
> = LazyLock::new(|| {
    Side::VALUES.map(|parent_side0| {
        Side::VALUES.map(|parent_side1| {
            let mut path_list: ArrayVec<[Side; 2], 2> = ArrayVec::new_with_default([Side::A; 2]);
            if parent_side0 == parent_side1 {
                return path_list;
            }
            for child_side0 in Side::iter() {
                if !child_side0.adjacent_to(parent_side0) || !child_side0.adjacent_to(parent_side1)
                {
                    // Child paths need to have both parts adjacent to parent paths.
                    continue;
                }
                for child_side1 in Side::iter() {
                    if child_side0 == child_side1 {
                        // Only look at child paths that obey shortlex rules
                        continue;
                    }
                    if child_side0.adjacent_to(child_side1)
                        && (child_side0 as usize) > (child_side1 as usize)
                    {
                        // Only look at child paths that obey shortlex rules
                        continue;
                    }
                    if !child_side1.adjacent_to(parent_side0)
                        || !child_side1.adjacent_to(parent_side1)
                    {
                        // Child paths need to have both parts adjacent to parent paths.
                        continue;
                    }
                    path_list.push([child_side0, child_side1]);
                }
            }
            path_list
        })
    })
});

trait GraphRef: AsRef<Graph> {
    fn length(&self, node: NodeId) -> u32;
    fn neighbor(&mut self, node: NodeId, side: Side) -> NodeId;
}

struct ImmutableGraphRef<'a> {
    graph: &'a Graph,
}

impl AsRef<Graph> for ImmutableGraphRef<'_> {
    #[inline]
    fn as_ref(&self) -> &Graph {
        self.graph
    }
}

impl GraphRef for ImmutableGraphRef<'_> {
    #[inline]
    fn length(&self, node: NodeId) -> u32 {
        self.graph.length(node)
    }

    #[inline]
    fn neighbor(&mut self, node: NodeId, side: Side) -> NodeId {
        self.graph.neighbor(node, side).unwrap()
    }
}

struct MutableGraphRef<'a> {
    graph: &'a mut Graph,
}

impl GraphRef for MutableGraphRef<'_> {
    #[inline]
    fn length(&self, node: NodeId) -> u32 {
        self.graph.length(node)
    }

    #[inline]
    fn neighbor(&mut self, node: NodeId, side: Side) -> NodeId {
        self.graph.ensure_neighbor(node, side)
    }
}

impl AsRef<Graph> for MutableGraphRef<'_> {
    #[inline]
    fn as_ref(&self) -> &Graph {
        self.graph
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peer_traverser_example() {
        let mut graph = Graph::new(1);
        //ensure_nearby(&mut graph, &Position::origin(), 6.0);
        let mut node = NodeId::ROOT;
        for side in [Side::B, Side::D, Side::C, Side::A] {
            node = graph.ensure_neighbor(node, side);
        }
        let mut traverser = PeerTraverser::new(node);
        while let Some(peer) = traverser.ensure_next(&mut graph) {
            println!(
                "Location: {:?}, {:?}, {:?}",
                graph.node_path(peer.node()),
                peer.path_from_base().collect::<Vec<_>>(),
                peer.path_from_peer().collect::<Vec<_>>()
            );
        }

        // TODO: Add a test that shows that path_from_base and path_from_peer
        // take you to the same node
    }
}
