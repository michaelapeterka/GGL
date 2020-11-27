rule[
 ruleID "r4: hv rearrangment_rev"
 left[
  edge[source 1 target 2 label "-"]
  edge[source 2 target 3 label "-"]
  edge[source 2 target 7 label "-"]
  edge[source 3 target 4 label "="]
  edge[source 3 target 8 label "-"]
  ]
 context[
  node[id 0 label "C"]
  node[id 1 label "C"]
  node[id 2 label "N"]
  node[id 3 label "C"]
  node[id 4 label "N"]
  node[id 5 label "*"]
  node[id 6 label "*"]
  node[id 7 label "H"]
  node[id 8 label "H"]
  edge[source 0 target 1 label "="]
  edge[source 0 target 4 label "-"]
  edge[source 0 target 5 label "-"]
 ]
 right[
  edge[source 1 target 3 label "-"]
  edge[source 3 target 2 label "#"]
  edge[source 4 target 7 label "-"]
  edge[source 4 target 8 label "-"]

 ]
]
