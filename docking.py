import urllib.request
import re
import os
import random
import math

thr = 4.
'''
def download_pdb(pdb_id):
    base_url = "https://files.rcsb.org/view/"
    file_name = pdb_id + ".pdb"
    target_url = base_url + file_name
    directory_path = "Desktop/xu/4750/csee012"
    pdb_file = directory_path + file_name
    if (file_name.lower() in os.listdir(directory_path)) or (
                file_name.upper() in os.listdir(directory_path)):
        print("\n" + "-----------------\n" +
              file_name + " is already in " + directory_path)  # +
    else:
        try:
            urllib.request.urlretrieve(target_url, pdb_file)
            print(file_name + " downloaded>>>>")
            #input("Pause01! Press enter to start again.")
        except(urllib.error.HTTPError or TimeoutError or urllib.error.URLError):
            print("ERROR: Network connection problem when "
                  "trying to download " + pdb_file + " from " + target_url + ". \nReturning to program...")
'''
def parse_pdb(pdb_id):
    file_name = pdb_id + ".pdb"
    directory_path = "C:/Structures/"
    pdb_file = directory_path + file_name
    atom = []
    het  = []
    na = 0 #no. of atom
    nh = 0 #no. of het atom
    #nr = 0 #no. of contact residues
    #caa = {} # contact amino acid
    #aap = {"ALA": 0,"ARG":0,"CYS":0,"ASP":0,"ASN":0,"GLU":0,"GLN":0,"MET":0,"LYS":0,"HIS":0,
    #       "SER": 0,"THR":0,"PHE":0,"TYR":0,"VAL":0,"GLY":0,"LEU":0,"ILE":0,"TRP":0,"PRO":0,"ALL":0,} # amino acid preferences
    with open(pdb_id+".pdb","r") as f:
        for line in f:
            if re.match('^ATOM.*', line):
                atom.append(" ")
                atom[na] = [" "]*7

                atom[na][0] = float(line[30:38].strip())  # x
                atom[na][1] = float(line[38:46].strip())  # y
                atom[na][2] = float(line[46:54].strip())  # z
                atom[na][3] = line[12:16]          # atom name
                atom[na][4] = line[21]             # chain ID
                atom[na][5] = line[17:20]          #res name
                atom[na][6] = int(line[22:26])          #res no.
                na += 1

            elif re.match('^HETATM.*STL', line):
                if line[21] == 'B':
                    continue
                het.append(" ")
                het[nh] = [" "] * 6
                het[nh][0] = float(line[30:38].strip()) # x
                het[nh][1] = float(line[38:46].strip()) # y
                het[nh][2] = float(line[46:54].strip()) # z
                het[nh][3] = line[12:16] #atom name
                het[nh][4] = line[21]
                nh += 1
            elif re.match('^(HEADER|COMPND|TITLE)',line):
                print(line)
                #input("pause01: press Enter to continue")
    #print("number of atoms:"+str(na)+" number of het atoms:"+str(nh))
    return atom, het

def initialization(c,m):
    #c = [[0.]*3 for i in range(4)] #list
    ix = random_n(m)
    iy = random_n(m)
    iz = random_n(m)
    for i in range(len(c)):
        c[i][0] += ix
        c[i][1] += iy
        c[i][2] += iz
    return c

def ratation(c,m):
    pi = math.pi
    a = random_n(m)/180.*pi #around z axis
    b = random_n(m)/180.*pi #around y axis
    g = random_n(m)/180.*pi #around x axis
    r = [[0]*3 for i in range(3)]
    r[0][0] = math.cos(a) * math.cos(b)
    r[0][1] = math.cos(a) * math.sin(b) * math.sin(g) - math.sin(a) * math.cos(g)
    r[0][2] = math.cos(a) * math.sin(b) * math.cos(g) + math.sin(a) * math.sin(g)
    r[1][0] = math.sin(a) * math.cos(b)
    r[1][1] = math.sin(a) * math.sin(b) * math.sin(g) + math.cos(a) * math.cos(g)
    r[1][2] = math.sin(a) * math.sin(b) * math.cos(g) - math.cos(a) * math.sin(g)
    r[2][0] = math.sin(b) * -1.
    r[2][1] = math.cos(b) * math.sin(g)
    r[2][2] = math.cos(b) * math.cos(g)


    n = [[0]*3 for i in range(len(c))]
    o = [[0]*3 for i in range(len(c))]
    ori = geometry_center(c)
    for i in range(len(c)):
        n[i][0] = c[i][0] - ori[0] #x'
        n[i][1] = c[i][1] - ori[1] #y'
        n[i][2] = c[i][2] - ori[2] #z'
    for i in range(len(c)):
        for j in range(3):
            o[i][j] = n[i][j]
    for i in range(len(c)):
        n[i][0] = r[0][0] * o[i][0] + r[0][1] * o[i][1] + r[0][2] * o[i][2]  # x"
        n[i][1] = r[1][0] * o[i][0] + r[1][1] * o[i][1] + r[1][2] * o[i][2]  # y"
        n[i][2] = r[2][0] * o[i][0] + r[2][1] * o[i][1] + r[2][2] * o[i][2]  # z"

    for i in range(len(c)):
        c[i][0] = n[i][0] + ori[0]  # x
        c[i][1] = n[i][1] + ori[1]  # y
        c[i][2] = n[i][2] + ori[2]  # z
    return c

def translation(c,m):
    dx = random_n(m) #delta x
    dy = random_n(m)
    dz = random_n(m)
    for i in range(len(c)):
        c[i][0] += dx
        c[i][1] += dy
        c[i][2] += dz
    return c

def geometry_center(c):
    v = [0.,0.,0.]
    for i in range(len(c)):
        v[0] += c[i][0] / float(len(c))  # x
        v[1] += c[i][1] / float(len(c))  # y
        v[2] += c[i][2] / float(len(c))  # z
    return v

def random_n(n):
    a = random.uniform(-1.0*n,n)
    return a

def box(a,c):
    (xmax,ymax,zmax) = (-1000.,-1000.,-1000.)
    (xmin,ymin,zmin) = (1000.,1000.,1000.)
    for i in range(len(a)):
        if a[i][0] > xmax:
            xmax = a[i][0]
        if a[i][1] > ymax:
            ymax = a[i][1]
        if a[i][2] > zmax:
            zmax = a[i][2]
        if a[i][0] < xmin:
            xmax = a[i][0]
        if a[i][1] > ymin:
            ymax = a[i][1]
        if a[i][2] > zmin:
            zmax = a[i][2]
    lx = xmax - xmin
    ly = ymax - ymin
    lz = zmax - zmin
    l = max(lx,ly,lz)
    return xmax,xmin,ymax,ymin,zmax,zmin
    #print(xmax,xmin,ymax,ymin,zmax,zmin,lx,ly,lz,l)
    #input("please enter to continue")

def centralization(a,xc,yc,zc):
    for i in range(len(a)):
        a[i][0] -= xc
        a[i][1] -= yc
        a[i][2] -= zc
    return a

def boundry(c,l):
    (xc, yc, zc) = geometry_center(c) #find the location of the center of drug
    dx = 0.
    dy = 0.
    dz = 0.
    if xc > 1/2. :
        dx = (1/2.-xc)
    if xc < 1/-2.:
        dx = (1/-2.-xc)
    if yc > 1/2. :
        dy = (1/2.-yc)
    if yc < 1/-2.:
        dy = (1/-2.-yc)
    if zc > 1/2. :
        dz = (1/2.-zc)
    if zc < 1/-2.:
        dz = (1/-2.-zc)
    for i in range(len(c)):
        c[i][0] += dx
        c[i][1] += dy
        c[i][2] += dz
    return c

def distance(v1,v2):
    d = ((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)**.5 #compute the Eular distance
    return d

def bounce(a,c):
    for i in range(len(c)):
        for j in range(len(a)):
            d = distance(a[j],c[i])
            if d < 1.:
                overlap = True
                break
    return overlap

def rec(c,old_c): #put c into old in order to record or recover
    for i in range(len(c)):
        old_c[i][0] = c[i][0]
        old_c[i][1] = c[i][1]
        old_c[i][2] = c[i][2]
    return old_c

def main():
    side_length = 80.
    fout = open("./2try.pdb",'w') #a
    (a, c) = parse_pdb("4qoh")
    #(xmax,xmin,ymax,ymin,zmax,zmin) = box(a,c) #to create a virtual box
    (xc, yc, zc) = geometry_center(a) #find the geometry center of protein
    a = centralization(a, xc, yc, zc) #move the protein to the origin
    (xc, yc, zc) = geometry_center(c)  # find the geometry center of compound
    c = centralization(c, xc, yc, zc) #move the drug to the drigin
    c = initialization(c,side_length/2.)  # move the drug to the origin
    for i in range(100): #30 steps
        for j in range(len(a)): #
            if a[j][4] != a[j-1][4] and j>0:
                fout.write("TER\n")
            line = "{0:6s}{1:5d} {2:4s} {3:3s} {4:1s}{5:4d}    {6:8.3f}{7:8.3f}{8:8.3f}\n".format('ATOM  ', j + 1,
                                                                                                  a[j][3], a[j][5],
                                                                                                  a[j][4], a[j][6], a[j][0],
                                                                                                  a[j][1], a[j][2])
            print(line)
            fout.write(line)
        fout.write("TER\n")
        for j in range(len(c)): #no. of atom
            line = "{0:6s}{1:5d} {2:4s} {3:3s} {4:1s}{5:4d}    {6:8.3f}{7:8.3f}{8:8.3f}\n".format('HETATM',j+1,c[j][3],'STL',c[j][4],1,c[j][0],c[j][1],c[j][2])
            print(line)
            fout.write(line)
        fout.write("TER\n")
        fout.write("END\n")
        overlap = True
        old_c = [[0.] * 6 for i in range(len(c))]
        rec(c,old_c) #
        while (overlap): #when protein and drug bounce together we will re do it
            c = translation(c,5.)
            c = ratation(c,30.)
            overlap = bounce(a,c)
            if overlap:
                c = old_c
        c = boundry(c,side_length)
        print(c)
    fout.close()
main()
