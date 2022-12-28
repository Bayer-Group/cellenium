-- Hierarchical functions for concepts
-- By convention, the functions don't include their own concept IDs in the result sets / paths,
-- only their children and parents.

drop function if exists concept_children_depth_limit;
create function concept_children_depth_limit(concept_inp concept, depth_limit int)
    returns setof concept
as
$$
with recursive concept_hier_cte(cid, parent_cid, path) as (select distinct ch.cid, null::integer, array []::int[]
                                                           from concept_hierarchy ch
                                                           where ch.cid = concept_inp.cid
                                                           union all
                                                           select ch.cid,
                                                                  ch.parent_cid,
                                                                  array_append(concept_hier_cte.path, ch.cid)
                                                           from concept_hierarchy ch,
                                                                concept_hier_cte
                                                           where ch.parent_cid = concept_hier_cte.cid
                                                             and coalesce(array_length(concept_hier_cte.path, 1), 0) < depth_limit)
select c.*
from concept_hier_cte
         join concept c on c.cid = concept_hier_cte.cid
where concept_hier_cte.cid != concept_inp.cid;
$$ language sql stable;


drop function if exists concept_all_children;
create function concept_all_children(concept_inp concept)
    returns setof concept
as
$$
select *
from concept_children_depth_limit(concept_inp, 10000)
$$ language sql stable;

drop function if exists concept_direct_children;
create function concept_direct_children(concept_inp concept)
    returns setof concept
as
$$
select *
from concept_children_depth_limit(concept_inp, 1)
$$ language sql stable;

drop function if exists concept_all_parents;
create function concept_all_parents(concept_inp concept)
    returns setof concept
as
$$
with recursive concept_hier_cte(cid, parent_cid) as (select distinct null::integer, ch.parent_cid
                                                     from concept_hierarchy ch
                                                     where ch.cid = concept_inp.cid
                                                     union all
                                                     select ch.cid, ch.parent_cid
                                                     from concept_hierarchy ch,
                                                          concept_hier_cte
                                                     where ch.cid = concept_hier_cte.parent_cid)
select distinct c.*
from concept_hier_cte
         join concept c on c.cid = concept_hier_cte.parent_cid;
$$ language sql stable;

drop function if exists concept_direct_parents;
create function concept_direct_parents(concept_inp concept)
    returns setof concept
as
$$
select distinct c.*
from concept_hierarchy ch
         join concept c on ch.parent_cid = c.cid
where ch.cid = concept_inp.cid;
$$ language sql stable;

drop function if exists concept_cid_array_to_codes;
CREATE FUNCTION concept_cid_array_to_codes(int[]) RETURNS text[] AS
$$
SELECT array_agg(concept.ont_code order by vals.i)
FROM unnest($1) with ordinality as vals(v, i)
         join concept on vals.v = concept.cid;
$$ LANGUAGE sql STABLE;


drop function if exists concept_all_parents_paths;
create or replace function concept_all_parents_paths(concept_inp concept)
    returns table
            (
                cid_path      int[],
                ont_code_path text[]
            )
as
$$
with all_paths_to_roots as (with recursive concept_hier_cte(cid, parent_cid, path) as (select distinct null::integer,
                                                                                                       concept_inp.cid parent_cid, -- (select cid from concept where ont_code ='58797008')
                                                                                                       array []::int[]
                                                                                       union all
                                                                                       select ch.cid,
                                                                                              ch.parent_cid,
                                                                                              array_prepend(ch.cid, concept_hier_cte.path)
                                                                                       from concept_hierarchy ch,
                                                                                            concept_hier_cte
                                                                                       where ch.cid = concept_hier_cte.parent_cid)
                            select distinct concept_hier_cte.path cid_path --, concept_cid_array_to_codes(concept_hier_cte.path) ont_code_path
                            from concept_hier_cte
                                     join concept c on c.cid = concept_hier_cte.cid)
select array_remove(return_the_path.cid_path, concept_inp.cid),
       concept_cid_array_to_codes(array_remove(return_the_path.cid_path, concept_inp.cid)) ont_code_path
from all_paths_to_roots return_the_path
where not exists(
    -- suppress return_the_path if there is another path that contains it, and has additional parents
        select 1
        from all_paths_to_roots compare_path
        where array_length(compare_path.cid_path, 1) > array_length(return_the_path.cid_path, 1)
          and compare_path.cid_path @> return_the_path.cid_path
    );
$$ LANGUAGE sql STABLE;


drop function if exists concept_test_potential_children;
create function concept_test_potential_children(concept_inp concept, potential_children_codes text[])
    returns setof concept
as
$$
select c.*
from concept_all_children(concept_inp) cac
         join concept c on cac.cid = c.cid and c.ont_code = any (potential_children_codes);
$$ LANGUAGE sql STABLE;

drop function if exists concept_children_paths_depth_limit;
create function concept_children_paths_depth_limit(concept_inp concept, depth_limit int)
    returns table
            (
                cid_path      int[],
                ont_code_path text[]
            )
as
$$
with recursive concept_hier_cte(cid, parent_cid, path) as (select distinct ch.cid,
                                                                           null::integer /*ch.parent_cid*/,
                                                                           array []::int[]
                                                           from concept_hierarchy ch
                                                           where ch.cid = concept_inp.cid -- (select cid from concept where ont_code ='58797008')
                                                           union all
                                                           select ch.cid,
                                                                  ch.parent_cid,
                                                                  array_append(concept_hier_cte.path, ch.cid)
                                                           from concept_hierarchy ch,
                                                                concept_hier_cte
                                                           where ch.parent_cid = concept_hier_cte.cid
                                                             and coalesce(array_length(concept_hier_cte.path, 1), 0) < depth_limit)
select distinct concept_hier_cte.path cid_path, concept_cid_array_to_codes(concept_hier_cte.path) ont_code_path
from concept_hier_cte
         join concept c on c.cid = concept_hier_cte.cid
where concept_hier_cte.cid != concept_inp.cid;
$$ LANGUAGE sql STABLE;


drop function if exists concept_all_children_paths;
create function concept_all_children_paths(concept_inp concept)
    returns table
            (
                cid_path      int[],
                ont_code_path text[]
            )
as
$$
select *
from concept_children_paths_depth_limit(concept_inp, 10000)
$$ language sql;



drop function if exists _semantic_order_impl;
drop function if exists _concept_hierarchy_minimum_trees_impl;
drop type if exists concept_weighted_parent cascade;
create type concept_weighted_parent as
(
    cid                          int,
    parent_cid                   int,
    semantic_relationship_weight real
);


/*
 This code follows "A new method to measure the semantic similarity of GO terms", Wang et al
 https://academic.oup.com/bioinformatics/article/23/10/1274/197095
 S, t, A, B are used according to the paper. The paper has additional complexity of considering the
 overlapping GO annotations of two genes, we don't have that here, we just calculate the semantic
 similarity of pairs of terms in an ontology.

 The terms are than clustered according to term distances (distance is the opposite of similarity),
 and the clustering yields a tree which leaves are in an order that has conceptual similar terms
 closer together.

 The ordered list of concepts can be used for e.g. a heatmap axis of cell types, diseases etc, keeping
 similar terms in proximity to each other.
 */

CREATE FUNCTION _semantic_order_impl(fullgraph concept_weighted_parent[], cids_to_order int[])
    RETURNS int[] AS
$$

import networkx as nx
from typing import List, Dict
from itertools import permutations
import numpy as np
import scipy
from scipy import spatial
from scipy import cluster


def calc_S_value(G, t, svalues: Dict):
    for a in G.pred[t]:
        # print(f"{t} - {a} {svalues[t]} {svalues[t]*0.8}  ({svalues.get(a,0)}) -> {max(svalues[t]*0.8, svalues.get(a,0))}")
        w = G.get_edge_data(a, t)
        w = w['weight']
        svalues[a] = max(svalues[t] * w, svalues.get(a, 0))
        calc_S_value(G, a, svalues)


def build_S_values_dict_for_concept(G, t) -> Dict:
    svalues = {t: 1.0}
    if t in G:
        calc_S_value(G, t, svalues)
    return svalues


def construct_graph(concept_weighted_parents: List):
    G = nx.DiGraph()
    for cp in concept_weighted_parents:
        G.add_nodes_from([cp['cid'], cp['parent_cid']])
        G.add_edge(cp['parent_cid'], cp['cid'], weight=cp['semantic_relationship_weight'])
    return G


def all_S_values_dicts(G, ts: List) -> Dict[int, Dict]:
    S_values_dicts = dict()
    for t in ts:
        S_values_dicts[t] = build_S_values_dict_for_concept(G, t)
    return S_values_dicts


def concept_distance(G, all_S_values_dicts, A, B):
    if A == B:
        return 0.0

    A_svalues = all_S_values_dicts[A]
    B_svalues = all_S_values_dicts[B]

    def sum_svalues(svalues: Dict, t: List):
        sum = 0.0
        for i in t:
            sum += svalues.get(i, 0.0)
        return sum

    common_terms = A_svalues.keys() & B_svalues.keys()
    all_terms = A_svalues.keys() | B_svalues.keys()
    common_terms_sum = sum_svalues(A_svalues, common_terms) + sum_svalues(B_svalues, common_terms)
    all_terms_sum = sum_svalues(A_svalues, all_terms) + sum_svalues(B_svalues, all_terms)
    # print(f"{A} - {B}  distance: {all_terms_sum/common_terms_sum}  (similarity: {common_terms_sum/all_terms_sum}))
    if common_terms_sum == 0.0:
        return 10000.0
    else:
        return all_terms_sum / common_terms_sum


def distance_matrix(G, S_values_dicts):
    # using plain 'S_values_dicts' works with the dict keys
    cid_to_matrixindex = {cid: i for i, cid in enumerate(S_values_dicts)}
    # print(cid_to_index)

    X = np.zeros(shape=(len(S_values_dicts), len(S_values_dicts)), dtype=float)
    for pair in permutations(S_values_dicts, 2):
        ai = cid_to_matrixindex[pair[0]]
        bi = cid_to_matrixindex[pair[1]]
        if ai < bi:
            X[ai, bi] = concept_distance(G, S_values_dicts, pair[0], pair[1])
            X[bi, ai] = X[ai, bi]
    return X, cid_to_matrixindex


def leaf_ordering_by_hierarchy_level(G, ts, Z, y):
    def calc_to_last_root(t):
        l = 0.0
        for ancestor in nx.ancestors(G, t):
            l = max(l, nx.shortest_path_length(G, ancestor, t, 'weight'))
        return l

    dist_to_root = {t: calc_to_last_root(t) for t in ts}
    Z_ordered = cluster.hierarchy.optimal_leaf_ordering(Z, y, metric=lambda A, B: -1 if dist_to_root[A] < dist_to_root[
        B] else 1)
    return Z_ordered


def hierarchical_clustering(G, cids_to_order, X):
    y = spatial.distance.squareform(X)
    Z = cluster.hierarchy.linkage(y, "complete")
    # can be plotted:
    # cluster.hierarchy.dendrogram(Z, orientation='left')
    # now swap the leaves/branches to have terms higher in the hierarchy also first in the dendrogram tree
    Z_ordered = leaf_ordering_by_hierarchy_level(G, cids_to_order, Z, y)
    # Z_ordered can also be plotted...
    rootnode = cluster.hierarchy.to_tree(Z_ordered)
    return rootnode


def ordered_leaves_cids(rootnode, cid_to_matrixindex) -> List[int]:
    cluster_ordered_matrix_indexes = rootnode.pre_order()
    matrixindex_to_cid = {i: cid for cid, i in cid_to_matrixindex.items()}
    ordered_cids = [matrixindex_to_cid[i] for i in cluster_ordered_matrix_indexes]
    return ordered_cids


def semantic_order(concept_weighted_parents: List[Dict], cids_to_order: List[int]):
    if len(cids_to_order) > 2000:
        raise Exception('semantic_order limited to up to 2000 terms')
    if len(cids_to_order) <= 1:
        return cids_to_order

    G = construct_graph(concept_weighted_parents)
    S_values_dicts = all_S_values_dicts(G, cids_to_order)
    X, cid_to_matrixindex = distance_matrix(G, S_values_dicts)
    rootnode = hierarchical_clustering(G, cids_to_order, X)
    # possible improvement:
    # swap the leaf and its siblings, depending on which has the higher level in ontology,
    # then go up to the higher-level branches and swap them as well, depending on which has the higher level nodes
    return ordered_leaves_cids(rootnode, cid_to_matrixindex)


return semantic_order(fullgraph, cids_to_order)

$$ LANGUAGE plpython3u
    IMMUTABLE
    SECURITY DEFINER
    PARALLEL SAFE;

/*
select * from _semantic_order_impl(array[
    (2,1,0.8),
    (3,2,0.8),
    (4,3,0.8),
    (5,3,0.8),
    (6,1,0.8)
    ]::concept_weighted_parent[], array[3,4,5,6]);
*/


drop function if exists concepts_in_semantic_order;
CREATE FUNCTION concepts_in_semantic_order(query_ontology text, ontology_codes text[]) RETURNS SETOF concept AS
$$
with query_concepts as (select distinct c.cid, c.ontid
                        from concept c
                                 join ontology o on c.ontid = o.ontid
                        where c.ont_code = any (ontology_codes)
                          and o.name = query_ontology),
     ordering_cids_and_all_parents as (with recursive concept_hier_cte(cid, parent_cid)
                                                          as (select distinct null::integer, query_concepts.cid parent_cid
                                                              from query_concepts
                                                              union all
                                                              select ch.cid, ch.parent_cid
                                                              from concept_hierarchy ch,
                                                                   concept_hier_cte
                                                              where ch.cid = concept_hier_cte.parent_cid)
                                       select distinct concept_hier_cte.cid, concept_hier_cte.parent_cid
                                       from concept_hier_cte
                                       where concept_hier_cte.cid is not null)
select c.*
from _semantic_order_impl(
             array(
                     select (ch.cid, ch.parent_cid, 0.8)::concept_weighted_parent
                     from ordering_cids_and_all_parents ch
                     union all
                     -- there may be ontology_codes which don't have any parents - add a fake no-op relationship
                     select (qc.cid, -1 * qc.cid, 0.1)::concept_weighted_parent
                     from query_concepts qc
                     where qc.cid not in (select cid from ordering_cids_and_all_parents)
                 ),
             (select array_agg(query_concepts.cid)
              from query_concepts)) semantic_order_result
         cross join unnest(semantic_order_result) with ordinality ordered_cid(cid, i)
         join concept c on c.cid = ordered_cid.cid
order by ordered_cid.i;
$$ LANGUAGE sql STABLE;


-- select ont_code, label from concepts_in_semantic_order('Cell Ontology', array ['CL_0000787','CL_0000791','CL_0000076','CL_0000694','CL_0002043','CL_0002044']);
-- select ont_code, label from concepts_in_semantic_order('Cell Ontology', array ['CL_0000791','CL_0002419', 'CL_0000084', 'CL_0000542']);
-- select ont_code, label from concepts_in_semantic_order('Cell Ontology', array ['CL_0000694','CL_2000019', 'CL_0000488', 'CL_0000210']);
-- select ont_code, label from concepts_in_semantic_order('Cell Ontology', array ['CL_0000791','CL_0002419', 'CL_0000084', 'CL_0000542','CL_0000694','CL_2000019', 'CL_0000488', 'CL_0000210']);
-- select ont_code, label from concepts_in_semantic_order('MeSH', array ['D006973', 'D003920', 'D002280']);


drop function if exists _concept_hierarchy_minimum_trees_impl, concept_hierarchy_minimum_trees, concept_hierarchy_minimum_trees_parents_lists;
drop type if exists concept_tree_element, concept_path, minimum_trees_result;

create type concept_tree_element as
(
    cid   int,
    level int
);
create type concept_path as
(
    cid         int,
    parent_cids int[]
);
create type minimum_trees_result as
(
    concept_tree_elements concept_tree_element[],
    concept_paths         concept_path[]
);

/*
 From the concepts given in cids_leaves (ts), find as few roots as possible to reach all given concepts.
 Then decide for the "most popular" path from root to concept in case multiple paths exist. Return the one
 path per concept. The structure can still consist of multiple trees, each representing a distinct set of
 the given concepts.
 */

drop function if exists _concept_hierarchy_minimum_trees_impl;
CREATE FUNCTION _concept_hierarchy_minimum_trees_impl(fullgraph concept_weighted_parent[], cids_leaves int[])
    RETURNS minimum_trees_result AS
$$

import networkx as nx
from typing import List, Dict
from itertools import permutations
import numpy as np
import pandas as pd


def construct_graph(concept_weighted_parents: List, cids_leaves: List):
    G = nx.DiGraph()
    for cp in concept_weighted_parents:
        G.add_nodes_from([cp['cid'], cp['parent_cid']])
        w = G.get_edge_data(cp['parent_cid'], cp['cid'], {'weight': 0})['weight']
        if w == 0:
            G.add_edge(cp['parent_cid'], cp['cid'], weight=-1)
        else:
            G[cp['parent_cid']][cp['cid']]['weight'] = w - 1
    # make sure all leaves are present in the graph
    G.add_nodes_from(cids_leaves)
    return G


def find_required_roots(G, ts):
    # https://stackoverflow.com/questions/62468287/finding-all-the-roots-inside-a-digraph-networkx
    all_roots = []
    for component in nx.weakly_connected_components(G):
        G_sub = G.subgraph(component)
        all_roots.extend([n for n, d in G_sub.in_degree() if d == 0])

    t_root_ts = []
    t_root_roots = []
    for t in ts:
        for r in all_roots:
            if nx.has_path(G, r, t):
                t_root_ts.append(t)
                t_root_roots.append(r)
            # try:
            #    # print(t, r, nx.shortest_path_length(G, r, t, weight='weight'))
            #    t_root_ts.append(t)
            #    t_root_roots.append(r)
            # except nx.NetworkXNoPath:
            #    pass
    df = pd.DataFrame({'t': t_root_ts, 'root': t_root_roots})
    # print(df.head())

    picked_roots = []
    while len(df) > 0:
        # print(df.head())
        frequent_root = df.root.mode().values[0]
        picked_roots.append(frequent_root)
        covered_t = df[df.root == frequent_root].t.values
        # find all remaining root-term relationships which are not covered yet
        df = df[~df.t.isin(covered_t)]
    return picked_roots


def build_trees(G, roots, ts):
    TG = nx.DiGraph()
    parent_paths = {}
    for t in ts:
        l = 0
        nodes = None
        for r in roots:
            if nx.has_path(G, r, t):
                li = nx.bellman_ford_path_length(G, r, t, weight='weight')
                if li < l:
                    nodes = nx.bellman_ford_path(G, r, t, weight='weight')
        if nodes:
            nx.add_path(TG, nodes)
            # same behavior as 'concept_all_parents_paths': start from concept but don't include it in the list, then work up to the root
            parent_paths[t] = list(reversed(nodes))[1:]
        else:
            TG.add_node(t)
            parent_paths[t] = list()
    return TG, parent_paths


def print_tree(TG, root):
    d = nx.dfs_successors(TG, root)
    concept_tree_elements = []

    def _level(n, level):
        # print(' ' * level + str(n))
        # plpy.info('.' * level + str(n))
        concept_tree_elements.append({'level': level, 'cid': n})
        for child in d.get(n, []):
            _level(child, level + 1)

    _level(root, 0)
    return concept_tree_elements


def minimum_trees(fullgraph, cids_leaves):
    G = construct_graph(fullgraph, cids_leaves)
    picked_roots = find_required_roots(G, cids_leaves)
    TG, parent_paths = build_trees(G, picked_roots, cids_leaves)

    concept_tree_elements = []
    for r in picked_roots:
        concept_tree_elements.extend(print_tree(TG, r))
    return concept_tree_elements, parent_paths


concept_tree_elements, parent_paths = minimum_trees(fullgraph, cids_leaves)

return {
    'concept_tree_elements': concept_tree_elements,
    'concept_paths': [
        {'cid': cid, 'parent_cids': parent_cids} for cid, parent_cids in parent_paths.items()
    ]
}

$$ LANGUAGE plpython3u
    IMMUTABLE
    SECURITY DEFINER
    PARALLEL SAFE;


drop function if exists concept_hierarchy_minimum_trees_parents_lists;
CREATE FUNCTION concept_hierarchy_minimum_trees_parents_lists(query_ontology text, ontology_codes text[])
    RETURNS
        -- minimum_trees_result
        --SETOF concept_tree_element

        TABLE
        (
            cid                  int,
            ont_code             text,
            label             text,
            parent_cids          int[],
            parent_ont_code_path text[]
        )
AS
$$
    begin
    if ontology_codes is not null then
    return query
with query_concepts as (select distinct c.cid, c.ontid
                        from concept c
                                 join ontology o on c.ontid = o.ontid
                        where c.ont_code = any (ontology_codes)
                          and o.name = query_ontology),
     ordering_cids_and_all_parents as (with recursive concept_hier_cte(cid, parent_cid)
                                                          as (select distinct null::integer, query_concepts.cid parent_cid
                                                              from query_concepts
                                                              union all
                                                              select ch.cid, ch.parent_cid
                                                              from concept_hierarchy ch,
                                                                   concept_hier_cte
                                                              where ch.cid = concept_hier_cte.parent_cid)
                                       select distinct concept_hier_cte.cid, concept_hier_cte.parent_cid
                                       from concept_hier_cte
                                       where concept_hier_cte.cid is not null)
select elems.cid, cid_c.ont_code, cid_c.label, elems.parent_cids, concept_cid_array_to_codes(elems.parent_cids)
from unnest((select concept_paths
             from
                 _concept_hierarchy_minimum_trees_impl(
                         array(
                                 select (ch.cid, ch.parent_cid, 0.8)::concept_weighted_parent
                                 from ordering_cids_and_all_parents ch
                             ),
                         (select array_agg(query_concepts.cid)
                          from query_concepts)
                     ))) as elems
         join concept cid_c on elems.cid = cid_c.cid;
    end if;
    end;
$$ LANGUAGE plpgsql STABLE;

/* similar to concept_all_parents_paths, but returns only one parent path (up to one root) per concept:
select * from concept_hierarchy_minimum_trees_parents_lists('NCIT', array ['C12391','C33645','C34125','C34192','C33105','C34168','C13063','C12403','C22600','C13056','C34167','C12377','C12423','C13300','C12313','C33797','C12372','C12382','C12434','C25726','C12445','C12400','C13049','C12666','C156591','C34170','C12801','C12464','C12416','C12393','C12433','C12679','C12388','C12366','C98275','C32425','C12470','C12745','C12387','C12710','C12405','C12419','C12410','C12439','C12390','C12376','C13363','C12401','C12468','C12756','C12414','C34169','C12727','C12669','C12432','C12381','C12392','C12736','C34172','C12379','C12782','C12431','C12311','C33177','C12263','C12469','C12935','C12412','C12802','C12386','C33209','C13272','C12811','C12272','C12472','C22673','C52975','C12389','C12971','C12404','C12428','C12426','C12415','C12623','C173496']);
 */

drop type if exists autocomplete_result;
create type autocomplete_result as (
ontology                       text,
                ont_code                     text,
                label                        text,
                label_highlight              text,
                is_synonym_of_preferred_term text
                                   );

drop function if exists autocomplete;
create function autocomplete(search_query text)
    returns setof autocomplete_result
    language sql stable
    as
    $$
with all_concepts_terms as (select cid, ontid, ont_code, label, label_tsvector, null::text is_synonym_of_preferred_term
                            from concept
                            union all
                            select cs.cid,
                                   c.ontid,
                                   c.ont_code,
                                   cs.synonym          "label",
                                   cs.synonym_tsvector "label_tsvector",
                                   c.label             is_synonym_of_preferred_term
                            from concept_synonym cs
                                     join concept c on cs.cid = c.cid),
     as_tsquery as (
         -- users may double-quote words to find them exactly, otherwise we assume a prefix search
         -- it is also possible to double-quote a phrase of multiple words
         select --to_tsquery('simple','cdk | cdk:*') q
                to_tsquery('simple', string_agg(
                        case
                            when right(split, 1) = '"' then replace(replace(split, '"', ''), ' ', ' <-> ')
                            else split || ' | ' || split || ':*' end,
                        ' & ')) q
                -- see https://stackoverflow.com/questions/4780728/regex-split-string-preserving-quotes
         from regexp_split_to_table(lower(regexp_replace(search_query, '[^\w|\-|\"]', ' ', 'g')),
                                    E'(?<=^[^"]*(?:"[^"]*\"[^"]*)*) (?=(?:[^"]*"[^"]*")*[^"]*$)') split),
     search_results as (select as_tsquery.q,
                               c.ontid,
                               c.ont_code,
                               ts_rank(c.label_tsvector, as_tsquery.q, 2) *
                               case
                                   when lower(c.label) =
                                        lower(regexp_replace(search_query, '[^\w|\-|\"]', ' ', 'g'))
                                       then 2 /*still need to boost the exact match...*/
                                   else 1
                                   end
                                   *
                               case
                                   when c.is_synonym_of_preferred_term is null
                                       then 2 /* synonym matches are ranked lower than preferred label matches (example KRAS/NRAS) */
                                   else 1
                                   end                            "rank",
                               c.label,
                               ts_headline(c.label, as_tsquery.q) label_highlight,
                               c.is_synonym_of_preferred_term
                        from all_concepts_terms c,
                             as_tsquery
                        where c.label_tsvector @@ as_tsquery.q),
     search_results_ordered as (select *,
                                       row_number() over (partition by ontid order by rank desc, label) ontology_row_number,
                                       row_number()
                                       over (partition by ontid, ont_code order by rank desc, label)    deduplicate_synonym_hits
                                from search_results)
select o.name ontology,
       sro.ont_code,
       sro.label,
       sro.label_highlight,
       sro.is_synonym_of_preferred_term
from search_results_ordered sro
join ontology o on sro.ontid = o.ontid
where sro.ontology_row_number <= 20 /* don't exceed n hits per domain */
  and sro.deduplicate_synonym_hits = 1
order by sro.rank desc, sro.ontology_row_number;
$$;



-- select * from autocomplete('bone');





