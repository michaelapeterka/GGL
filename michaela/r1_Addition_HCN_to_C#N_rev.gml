rule[
 ruleID "r1: Addition HCN to C#N_rev"
 left[
  node[id 1 label "C"]
  node[id 4 label "H"]
  edge[source 1 target 2 label "-"]
  edge[source 2 target 3 label "="]
  edge[source 3 target 4 label "-"]
 ]
 context[
  node[id 0 label "N"]
  node[id 2 label "C"]
  node[id 3 label "N"]
  node[id 5 label "*"]
  edge[source 0 target 1 label "#"]
  edge[source 2 target 5 label "-"]
 ]
 right[
  node[id 1 label "C-"]
  node[id 4 label "H+"]
   edge[source 2 target 3 label "#"]
 ]
]
