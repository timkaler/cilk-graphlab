import random
import sys

edge_dict = dict()

def add_edge(i,j,weight):
  if (i,j) not in edge_dict:
    edge_dict[(i,j)] = True
    #print str(i)+"\t"+str(j)+"\t"+str(weight)+"\n"
    print str(i)+"\t"+str(j)+"\n"
    return True
  else:
    return False

def main():
  random.seed(53)
  if len(sys.argv) != 3:
    print "Usage: graph_gen.py <vertices> <degree>"
    sys.exit(1)
  V = int(sys.argv[1])
  degree = int(sys.argv[2])
  for i in range(0, V):
    n = dict()
    for j in range(0, degree/2):
      n[random.randint(0,V-1)] = True
    for j in n.keys():
      add_edge(i,j, random.random())

if __name__ == "__main__":
    main()
