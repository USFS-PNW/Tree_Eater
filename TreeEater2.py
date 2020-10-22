#!/usr/bin/env python3.7
import io
import cmath
import numpy
import os
import sys
from numbers import Number
os.chdir(os.path.dirname(sys.argv[0]))

# ----- USER CONTROLS -----
# mcbark and mcwood set the assumed moisture content for bark and wood respectively.
# See Neo-Processor documentation for full description of moisture content settings.
# Set target to list of cut tree records.
mcbark = 0
mcwood = 0
rx = "test"
tar = "MerchTrees.txt"

# ----- END OF USER CONTROLS -----

def gendib(tree, height, mdib):
    """gendib takes a commercial tree entry and determines the diameter inside bark at a given height.
    PICO equation taken from the corrected Garber and Maguire 2002, all others from Hann 2016."""
    if tree[9] != 5:
        cs = tree[9]
        pdib = (mdib[cs][0] * (tree[1] ** mdib[cs][1])) * cmath.exp((mdib[cs][2]) * ((1 - tree[4]) ** 0.5))
        hd = (tree[2] - 4.5) / tree[1]
        height = height - 4.5
        rh = height / (tree[2] - 4.5)
        wlt = (mdib[cs][7] * tree[11] - 4.5) / (tree[2] - 4.5)
        i1 = 0
        i2 = 0
        if rh >= 0 and wlt >= rh:
            i1 = 0
        if wlt < rh and rh <= 1:
            i1 = 1
        if wlt <= 0:
            i2 = 0
        if wlt > 0:
            i2 = 1
        jp1 = (rh - 1) / (wlt - 1)
        jp2 = (wlt - rh) / (wlt - 1)
        x1 = mdib[cs][3] + (mdib[cs][4] * cmath.exp(mdib[cs][5] * (hd ** 2)))
        x2 = mdib[cs][6]
        z0 = 1.0 - rh + i2 * (rh + i1 * (jp1 * (1.0 + jp2) - 1.0 )) - (rh - 1) * (rh - i2 * rh)
        z1 = (i2 * (rh + i1 * (jp1 * (rh + (wlt * jp2)) - rh)) - ((rh - 1) * (rh - i2 * rh))) 
        z2 = i2 * ((rh**2) + i1 * (jp1 * wlt * (2 * rh - wlt + (wlt * jp2)) - rh**2))
        dib = pdib * (z0 + (x1 * z1) + (x2 * z2))
        dib = dib.real
    if tree[9] == 5:
        cs = tree[9]
        mdbh = tree[1] * 2.54
        mht = tree[2] / 3.28084
        height = height / 3.28084
        z = height / mht
        q = 1 - z**0.5
        p = 1.37 / mht
        x = q / (1 - p**0.5)
        d = (mdbh**2) / mht
        c = (mdib[cs][1] * cmath.asin(q)) + (mdib[cs][2] * cmath.log(x)) + (mdib[cs][3] * cmath.log(z)) + (mdib[cs][4] * cmath.exp(z)) + (mdib[cs][5] * d) + (mdib[cs][6] * cmath.exp(z) * d)
        dib = mdib[cs][0] * mdbh * (x**c)
        dib = dib.real
        dib = dib / 2.54
    return dib

def dibtree(tree, mdib):
    """Takes a tree and populates a list with dibs. First nine indexes (0-8) should be empty.
    Every other index will be filled with the appropriate dib for that height until dib < 4 or
    until it reaches the actual height of the tree"""
    count = 9
    diblist = [0,0,0,0,0,0,0,0,0]
    while count < tree[3]:
        dib = gendib(tree,count,mdib)
        if dib >= 4:
            diblist.append(dib)
            count = count + 1
        if dib < 4:
            count = 999
    return diblist

def pricelog(topdib, length, species, slist, clist):
    """Takes a log with a given topdib, length, and species and returns a price for that log."""
    bf = slist[round(length/4)][round(topdib)]
    mbf = bf / 1000
    if topdib >= 6 and topdib <= 8:
        pbracket = 8
    if topdib > 8 and topdib <= 14:
        pbracket = 14
    if topdib > 14 and topdib <= 22:
        pbracket = 22
    if topdib > 22 and topdib <= 24:
        pbracket = 24
    p = mbf * clist[species][pbracket]
    return p

def gennodes(maxht):
    """Generates an empty node list for use in bucking."""
    count = 0
    nodes = list()
    while count <= maxht + 1:
        nodes.append(0)
        count = count + 1
    return nodes

def evalnode(tree, diblist, nodes, cnode, llist, maxht):
    """Evaluates all potential logs coming from a single node and saves the value of those logs in the forward node.
    Carries values current node forward where appropriate."""
    for length in llist:
        if cnode + length <= maxht:
            tnode = cnode + length + 1
            if diblist[tnode-1] >= 6:
                tval = pricelog(diblist[tnode-1], length, tree[9], slist, clist)
                pathval = tval + nodes[cnode]
            else:
                pathval = nodes[cnode]
            if nodes[tnode] < pathval:
                nodes[tnode] = pathval
    return nodes

def evaltree(tree, mdib, clist, slist, llist):
    """Takes a tree and generates a dibtree, then bucks that tree and values the outputs.
    Returns a nodelist with the optimal value for in each node.
    Whole tree value in $/cf is stored in index 0 of the nodelist."""
    diblist = dibtree(tree, mdib)
    maxht = len(diblist) - 1
    maxnode = maxht - 8
    nodes = gennodes(maxht)
    count = 1
    while count <= maxnode:
        nodes = evalnode(tree, diblist, nodes, count, llist, maxht)
        count = count + 1
    bestval = max(nodes)
    nodes[0] = bestval * (1 - tree[5])
    nodes[1] = bestval
    return nodes

def evalpulp(tree, nodes, mdib):
    """Takes a bucked tree and finds the lowest node with the optimal value.
    Uses that node to calculate the remaining volume to a four inch top and (if possible) create a pulp log.
    Pulp logs must be at least eight feet long and may go down to a diameter inside bark of four inches.
    Returns the cubic foot volume of any pulp log found in [0], non-log pulp in [1]."""
    optval = nodes[1]
    maxht = len(nodes) - 1
    count = 9
    pulp = [0,0]
    pulplen = 0
    while count <= maxht:
        if nodes[count] != optval:
            count = count + 1
        elif nodes[count] == optval:
            optht = count
            count = maxht + 1
            pulplen = maxht - optht
    if pulplen >= 8:
        pulptop = gendib(tree,maxht,mdib)
        pulpbot = gendib(tree,optht,mdib)
        r1 = pulptop / 24
        r2 = pulpbot / 24
        pulp[0] = (numpy.pi * pulplen * ((r1**2) + (r2**2) + (r1*r2))) / 3
        pulp[0] = pulp[0] * (1 - tree[5])
    if pulplen < 8:
        pulptop = gendib(tree,maxht,mdib)
        pulpbot = gendib(tree,optht,mdib)
        r1 = pulptop / 24
        r2 = pulpbot / 24
        pulp[1] = (numpy.pi * pulplen * ((r1**2) + (r2**2) + (r1*r2))) / 3
        pulp[1] = pulp[1] * (1 - tree[5])
    return pulp            

def calcbarkwt(volume, spcd, mc, gwp):
    """Takes a volume and species code and returns the green weight of bark in pounds. Requires gwp.
    Positive moisture content values set moisture content to that percent.
    Negative moisture content values will reduce FIADB moisture content values by that value (to a minimum of 10%).
    Set moisturecontent to zero to use FIADB moisture content values.
    Will print an error and return a weight of zero if species is not present in gwp."""
    barkwt = 0
    if spcd not in gwp:
        print("Species code " + str(spcd) + " not found in green weight parameters.")
    if spcd in gwp and mc > 0:
        mc = float(mc)
        barkwt = volume * 62.4 * gwp[spcd][2] * (gwp[spcd][0]/(100 + gwp[spcd][0])) * (1 + mc/100)
    if spcd in gwp and mc == 0:
        barkwt = volume * 62.4 * gwp[spcd][2] * (gwp[spcd][0]/(100 + gwp[spcd][0])) * (1 + gwp[spcd][4]/100)
    if spcd in gwp and mc < 0:
        mc = float(mc)
        mc = gwp[spcd][4] + mc
        if mc >= 10.0:
            barkwt = volume * 62.4 * gwp[spcd][2] * (gwp[spcd][0]/(100 + gwp[spcd][0])) * (1 + mc/100)
        else:
            barkwt = volume * 62.4 * gwp[spcd][2] * (gwp[spcd][0]/(100 + gwp[spcd][0])) * 1.1
    return barkwt

def calcwoodwt(volume, spcd, mc, gwp):
    """Takes a volume and species code and returns the green weight of wood in pounds. Requires gwp.
    Positive moisture content values set moisture content to that percent.
    Negative moisture content values will reduce FIADB moisture content values by that value (to a minimum of 10%).
    Set moisturecontent to zero to use FIADB moisture content values.
    Will print an error and return a weight of zero if species is not present in gwp."""
    woodwt = 0
    if spcd not in gwp:
        print("Species code " + str(spcd) + " not found in green weight parameters.")
    if spcd in gwp and mc != 0:
        mc = float(mc)
        woodwt = volume * 62.4 * gwp[spcd][1] * (1 - gwp[spcd][0]/(100 + gwp[spcd][0])) * (1 + mc/100)
    if spcd in gwp and mc == 0:
        woodwt = volume * 62.4 * gwp[spcd][1] * (1 - gwp[spcd][0]/(100 + gwp[spcd][0])) * (1 + gwp[spcd][3]/100)
    if spcd in gwp and mc < 0:
        mc = float(mc)
        mc = gwp[spcd][3] + mc
        if mc >= 10.0:
            woodwt = volume * 62.4 * gwp[spcd][1] * (1 - gwp[spcd][0]/(100 + gwp[spcd][0])) * (1 + mc/100)
        else:
            woodwt = volume * 62.4 * gwp[spcd][1] * (1 - gwp[spcd][0]/(100 + gwp[spcd][0])) * 1.1
    return woodwt

def incut(rx, woodmc, barkmc, filename, eqdict, gwp, mdib, clist, slist, llist):
    """incut imports a cutlist file and creates a pair of dictionaries, one of bucked tree records [0],
    one of per acre outputs by stand [1]. Requires unpacking."""
    trx = dict()
    orx = dict()
    tcount = 0
    scount = 0
    fhandle = io.open(filename, "r")
    for line in fhandle:
# Initial input processing and cut tree creation
        tlist = line.strip().split("\t")
        a = int(tlist[0])
        b = float(tlist[1])
        c = float(tlist[2])
        d = float(tlist[3])
        if d < 1:
            d = c
        e = float(tlist[4]) * 0.01
        f = float(tlist[5]) * 0.01
        g = float(tlist[6])
        h = float(tlist[7])
        i = str(tlist[8])
        if a not in eqdict:
            j = 0
        if a in eqdict:
            j = eqdict[a]
        if isinstance(tlist[9],Number):
            k = int(tlist[9])
        else:
            k = str(tlist[9])
        cht = d * e
        l = d - cht
        m = int(tlist[10])
        n = rx
        tree = [a, b, c, d, e, f, g, h, i, j, k, l, m, n]
# Cut tree testing (species, size, and ability to create a merch log) and creation of Tprice, MerchVol, and PulpVol
        eval = 1
        if tree[9] == 0:
            eval = 0
        if tree[1] < 6 or tree[1] > 24:
            eval = 0
        if eval == 1:
            f9 = gendib(tree, 9.0, mdib)
            if f9 < 6:
                eval = 0
        if eval == 0:
            tree.extend([0,0,tree[7],0])
        if eval == 1:
            nodes = evaltree(tree, mdib, clist, slist, llist)
            pulp = evalpulp(tree, nodes, mdib)
            tree.append(nodes[0])
            tree.append(pulp[1])
            tree.append(pulp[0])
            tree.append(((tree[7] * (1 - tree[5])) - pulp[1]) - pulp[0])
        if tree[16] > 0:
            pbark = calcbarkwt(tree[16],tree[0],barkmc,gwp)
            pwood = calcwoodwt(tree[16],tree[0],woodmc,gwp)
            ptotal = pbark + pwood
            ptons = ptotal / 2000
            tree.append(ptons)
        else:
            tree.append(0.0)
        if tree[17] > 0:
            mbark = calcbarkwt(tree[17],tree[0],barkmc,gwp)
            mwood = calcwoodwt(tree[17],tree[0],woodmc,gwp)
            mtotal = mbark + mwood
            mtons = mtotal / 2000
            tree.append(mtons)
        else:
            tree.append(0.0)      
# Storing cut tree in appropriate stand dictionary and creating that dictionary if it does not already exist
        if tree[8] not in trx:
            trx[tree[8]] = dict()
            scount = scount + 1
        trx[tree[8]][tree[10]] = tree
        tcount = tcount + 1
# Adds cut tree values to existing Per Acre Outputs,
        if tree[8] not in orx:
            orx[tree[8]] = [rx,tree[8],0,0,0,0,0]
        orx[tree[8]][2] = orx[tree[8]][2] + (tree[14] * tree[6])
        orx[tree[8]][3] = orx[tree[8]][3] + (tree[17] * tree[6])
        orx[tree[8]][4] = orx[tree[8]][4] + (tree[16] * tree[6])
        orx[tree[8]][5] = orx[tree[8]][5] + (tree[19] * tree[6])
        orx[tree[8]][6] = orx[tree[8]][6] + (tree[18] * tree[6])
    fhandle.close()
    print(str(tcount) + " cut tree records across " + str(scount) + " stands loaded for Rx " + str(rx))
    return [trx, orx]

def printout(printlist, fhandle, delimiter):
    """Writes a list out to an external file using delimited fields. Does
    not open or close the file handle. Does not return an object."""
    l = len(printlist)
    c = 1
    for index in printlist:
        out = str(index)
        fhandle.write(out)
        if c < l:
            fhandle.write(delimiter)
        c = c + 1
    fhandle.write("\n")
    return

# ----- DICTIONARY DEFINITIONS AND LOADING -----

# "s[x]" dictionaries store board foot volume of a log with a given top dib and length
# se and s1 should always be empty
# Called with slist[length][topdib]
# Length is in four foot segments

s1 = dict()         
s2 = dict()
s3 = dict()
s4 = dict()
s5 = dict()
s6 = dict()
s7 = dict()

slist = [0, s1, s2, s3, s4, s5, s6, s7]
llist = [8,12,16,20]

fhandle = io.open("scribtable.txt", "r")
next(fhandle)   #skip first line with headers
for line in fhandle:
    tlist = line.strip().split("\t")
    k = int(tlist[0])
    a = float(tlist[1])
    b = float(tlist[2])
    c = float(tlist[3])
    d = float(tlist[4])
    e = float(tlist[5])
    f = float(tlist[6])
    s2[k] = a
    s3[k] = b
    s4[k] = c
    s5[k] = d
    s6[k] = e
    s7[k] = f
fhandle.close()

# "c[x]" dictionaries hold prices per mbf by species and top dib
# ce should always be empty
# Called with clist[tree eater species code][dclass]
# PSME = 1
# CADE = 2
# PIPO = 3
# PILA = 4
# PICO = 5
# TFIR = 6

c1 = {8:481.25, 14:513.75, 22:526.25, 24:532.50}
c2 = {8:635.00, 14:635.00, 22:635.00, 24:635.00}
c3 = {8:282.50, 14:316.25, 22:345.00, 24:377.50}
c4 = {8:275.00, 14:295.00, 22:330.00, 24:345.00}
c5 = {8:342.50, 14:342.50, 22:342.50, 24:342.50}
c6 = {8:377.50, 14:391.25, 22:410.00, 24:416.25}

clist = [0, c1, c2, c3, c4, c5, c6]

# eqict takes FIA species codes and returns the tree eater species code
eqdict = {}
fhandle = io.open("SpeciesGroups.txt", "r")
next(fhandle)   #skip first line with headers
for line in fhandle:
    (key, val) = line.split()
    eqdict[int(key)] = int(val)
fhandle.close()

# dib lists contain the parameters neccessary for dib calculations
# See preptable.xlsx for list of indexes
# Call with mdib[species code][index]

dib1 = list()         
dib2 = list()
dib3 = list()
dib4 = list()
dib5 = list()
dib6 = list()

mdib = [0, dib1, dib2, dib3, dib4, dib5, dib6]

fhandle = io.open("dibparms.txt", "r")
parc = 1
for line in fhandle:
    tlist = line.strip().split("\t")
    mdib[parc].append(float(tlist[0]))
    mdib[parc].append(float(tlist[1]))
    mdib[parc].append(float(tlist[2]))
    mdib[parc].append(float(tlist[3]))
    mdib[parc].append(float(tlist[4]))
    mdib[parc].append(float(tlist[5]))
    mdib[parc].append(float(tlist[6]))
    mdib[parc].append(float(tlist[7]))
    parc = parc + 1
fhandle.close()

# Loading parameters for calculating green weight values

gwp = dict()
fhandle = io.open("gwtparms.txt", "r")
next(fhandle)   #skip first line with headers
for line in fhandle:
    tlist = line.strip().split("\t")
    a = int(tlist[0])
    b = float(tlist[1])
    c = float(tlist[2])
    d = float(tlist[3])
    e = float(tlist[4])
    f = float(tlist[5])    
    gwp[a] = [b, c, d, e, f]
fhandle.close()

cut = incut(rx,mcbark,mcwood,tar,eqdict,gwp,mdib,clist,slist,llist)
trx = cut[0]
orx = cut[1]
t_lab = "priced_t.txt"
s_lab = "priced_s.txt"

fhandle = open(t_lab,"w")
# Print headers before printing trees
t_headers = ["SPCD", "DBH", "HEIGHT", "TRUN_HT", "PCT_CR", "MDEFECT", "TPA", "MERCH_VOL", "STAND_ID", "TESC", \
             "TREE_ID", "HCB", "YEAR_CUT", "RX", "T_PRICE", "S_PULP_CF", "LOG_PULP_CF", "SAW_VOL_CF", "LOG_PULP_GT", \
             "SAW_WT_GT"]
t = 0
delimiter = "\t"
for h in t_headers:
    fhandle.write(h)
    t = t + 1
    if t < len(t_headers):
        fhandle.write(delimiter)
fhandle.write("\n")
count = 0
for stand in trx:
    for tree in trx[stand]:
        printout(trx[stand][tree],fhandle,"\t")
        count = count + 1
fhandle.close()
print(str(count) + " tree records written to " + t_lab)

fhandle = open(s_lab,"w")
count = 0
for stand in orx:
    printout(orx[stand],fhandle,"\t")
    count = count + 1
fhandle.close()
print(str(count) + " stand records written to " + s_lab)

print("DONE.")


    
