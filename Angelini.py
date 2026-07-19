import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
np.seterr(divide='ignore', invalid='ignore')


'''
------------------------------------------------------------------------------------------------------------------------
                                             Algorithm BFS-based and DFS-based implementation
------------------------------------------------------------------------------------------------------------------------
 @author: Evangelidakis Leandros
 @School of Applied Mathematics and Physical Sciences 

 For the paper tha accompanies the algorithm of Angelini et al. see follow this link:
 https://www.emis.de/journals/JGAA/accepted/2012/Angelini+2012.16.1.pdf

 -----------------------------------------------------------------------------------------------------------------------
 The graph internal data structures are based on an adjacency list representation and implemented using Python dictionary data structures. 
 The graph adjacency structure is implemented as a Python dictionary of dictionaries; the outer dictionary is keyed by nodes to values that 
 are themselves dictionaries keyed by neighboring node to the edge attributes associated with that edge.
------------------------------------------------------------------------------------------------------------------------
'''

def getGridArea(tree,root_node='A',display=False,save=False,algo="bfs",view=None,filename="tree1.png",w_labels=False,n_size=25,n_color='black',l_color='black'):
  
  g = tree.copy()
  graph = dict()

# The adjacency list of the input tree is stored in a dictionary of the form <node> : <list>
  for line in nx.generate_adjlist(g):
    ln = line.split(" ")
    graph.update({ ln[0]: ln[1:] })

#  map each node with its subtrees: node |-> T1(node),...
  subtree_map = dict()
  newg = g.copy()
  for k in list(g.nodes()):
    newg = g.copy()
    newg.remove_node(k)
    if graph[k] == [] :
      subtree_map.update({k: [k]})
    else:
      subtree_map[k] = []
      for node in graph[k]:
        v = list(nx.bfs_tree(newg, node).nodes())
        subtree_map[k].append(v)
        newg.remove_node(node)

#------------------------------------------------------------------------------------------------------------------------
#   Calculate the elements of the Stern-Brocot tree:
#
#   Stern-Brocot-Tree elements are of the form of a(n)/a(n+1), where
#   a(n) is derived from the Stern-Brocot Sequence as described in the
#   Online Encyclopedia of Integer Sequences (OEIS), sequence label: A002487
#------------------------------------------------------------------------------------------------------------------------
  n = len(graph.keys())

  if algo=="bfs":
    v = np.zeros(2*n)
    v[0] = 1
    for i in range(n):
      v[2 * i] = v[i]
      v[2 * i + 1] = v[i] + v[i + 1]

    x = v[1:n]
    y = v[2:n+1]

  elif algo=="dfs":
    y = np.zeros(n-1)
    x = np.ones(n-1)
    for i in range(n-1):
        y[i] = i+1

#------------------------------------------------------------------------------------------------------------------------    
# -Calculate the ratios y/x and order them by increasing value.
# -Combine them into an iterable of tuples object with zip() function.
# -Map each element of the ratios array with each tuple of the list,
#  in a dictionary called Sorted_Map
#------------------------------------------------------------------------------------------------------------------------  

  d = dict(zip(list(np.divide(y, x)),list(zip(y,x)))) 
  Sorted_Map = dict(sorted(d.items()))
  sorted_ratios = list(Sorted_Map.keys())        

#------------------------------------------------------------------------------------------------------------------------
#    -Consider the subtrees of root node and assign to each subtree
#     their corresponding sequences from the sorted S-B sequence array
#    -Also map for each roots adjacent node, his corresponding coordinates (y,x)
#------------------------------------------------------------------------------------------------------------------------


  subtrees_seq = dict() # maps each node subtrees with a corresponding sequence from the sorted S-B ratios
  coordmap = dict() # maps each vertex with its (y,x) coordinates (not relative to parent)
  coordmap.update({root_node: [0, 0]})
  T_u = dict()  # maps each vertex u with its tree's T(u)- sequence
  T_u.update({root_node:list(sorted_ratios)})
  for node in graph.keys():       # for every node
    if graph[node] != []:  
      #print(node)
      Tlength = []
      Tlength.append(0)
      for j in range(len(subtree_map[node])): # For each subtree of this node, fill the Tlength array containing their lengths
        Tlength.append(len(list(subtree_map[node])[j]))
      for j in range(len(subtree_map[node])): # For each subtree, calculate the start and end indexes
        start = 1 + sum(Tlength[:j + 1])
        end = sum(Tlength[:j + 2])
        subtrees_seq.update({''.join(subtree_map[node][j]): T_u[node][start - 1:end]}) # Map the subtree with a sequence
      for child in graph[node]: # For every child of the current node
        for k in subtree_map[node]:
          if child in k:
            coordmap.update({child:Sorted_Map[list(subtrees_seq[''.join(k)])[-1]]})  # Map the child with his coords
            T_u.update({child:subtrees_seq[''.join(k)]}) 
                           # Map the child with his corresponding sequence


  #------------------------------------------------------------------------------------------------------------------------
  # dictionary : grid_map
  # Points each node with its (y,x)-coordinates relative to parent node
  #------------------------------------------------------------------------------------------------------------------------
  

  grid_map = dict()
  grid_map.update({root_node : [0,0] })
  gvertices = list(graph.keys())
  for parent in gvertices:
    for child in graph[parent]:
      if child not in grid_map.keys():
        xparent = list(grid_map[parent])[1]
        yparent = list(grid_map[parent])[0]
        x = list(coordmap[child])[1]
        y = list(coordmap[child])[0]
        grid_map.update({child : [y+yparent,x+xparent]})
      else:
        continue

 # Make a dictionary of the form: { node : (x,y) coordinates }

  gridpos = dict()
  for node in g.nodes():
    #gridpos.update({node: np.array([grid_map[node][1], grid_map[node][0]]).astype(int)})
    gridpos.update({ node : (grid_map[node][1], grid_map[node][0]) })


  data = pd.DataFrame.from_dict(gridpos).to_numpy()
  x = data[0]
  y = data[1]
  maxx = int(max(x))+1
  maxy = int(max(y))+1
  gridArea = n**2
  graphArea = maxx*maxy


  #------------------------------------------------------------------------------------------------------------------------
  #                                           Plot Graph usign matplotlib
  #------------------------------------------------------------------------------------------------------------------------
  
  if (display == True) or (save == True):
    plt.figure(figsize=(7,7))
    nx.draw_networkx(g, pos=gridpos, with_labels=w_labels, node_size = n_size,arrows=False,node_color=n_color,font_color=l_color)
    plt.grid(color="gray")
    #plt.title('Algorithm '+algo.upper()+"-based",fontsize=14,fontweight='bold')
    if algo == 'bfs':
      plt.title('Angelini et al. BFS-based\n\nGrid Size: '+str(maxx)+' x '+str(maxy)+' ('+str(n)+' nodes)') 
    if algo == 'dfs':
      plt.title('Angelini et al. DFS-based\n\nGrid Size: '+str(maxx)+' x '+str(maxy)+' ('+str(n)+' nodes)') 

    ax = plt.gca()
    ax.set_ylim([-1,maxy+1])
    ax.set_xlim([-1,maxx+1])

    plt.xticks(np.arange(0,maxx+1,1)) 
    plt.yticks(np.arange(0,maxy+1,1)) 
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    #ax.set_aspect(1)
    #plt.tight_layout(pad=3)
    if save == True:
      plt.savefig(filename,bbox_inches='tight')
    if display == True:
      plt.show()
    
    plt.close()


  return [graphArea,maxx,maxy,gridpos] 
  