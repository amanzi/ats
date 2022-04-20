import sys,os

os.rename("Mesh.txt", "Mesh_with_spaces.txt")

fout = open("Mesh.txt",'w')
with open("Mesh_with_spaces.txt",'r') as fin:
    for line in fin:
        fout.write("\t".join(line.strip().split())+"\n")

fout.close()
