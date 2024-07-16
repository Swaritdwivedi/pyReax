def read_ffield(file_ffield):
    with open(file_ffield,"r") as f_ffield:
        gen_type = ['boc1','boc2','coa2','trip4','trip3','kc2','ovun6','trip2',
                 'ovun7','ovun8','trip1','swa','swb','n.u.','val6','lp1',
                 'val9','val10','n.u.','pen2','pen3','pen4','n.u.','tor2',
                 'tor3','tor4','n.u.','cot2','vdw1','cutoff','coa4','ovun4',
                 'ovun3','val8','acut','hbtol','n.u.','n.u.','coa3']
        numLine=1
        para={}
        for line in f_ffield:
            if numLine==1:
                remark=line
            elif "!" in line and ";" not in line:
                #print(" ".join(line.split("!")[1].split("\n")[0].split()))
                para[" ".join(line.split("!")[1].split("\n")[0].split())]=float(line.split("!")[0].split()[0])
                #print(line.split("!")[1])
            elif "!" in line and ";" in line:
                if "." not in line:
                    para[" ".join(line.split("!")[1].split(";")[0].split())]=float(line.split("!")[0].split()[0])
                else:
                    para[" ".join(line.split("!")[1].split(";")[0].split(".")[0].split())]=float(line.split("!")[0].split()[0])
            else:
                pass

            numLine+=1
        nGen=para["Number of general parameters"]
        nAtom=para["Nr of atoms"]
        nBond=para["Nr of bonds"]
        nOff=para["Nr of off-diagonal terms"]
        nAng=para["Nr of angles"]
        nTor=para["Nr of torsions"]
        nHB=para["Nr of hydrogen bonds"]
        numLine_gen=[int(3),int(3+nGen-1)]
        numLine_atom=[int(numLine_gen[1]+5),int(numLine_gen[1]+5+nAtom*4-1)]
        numLine_bond=[int(numLine_atom[1]+3),int(numLine_atom[1]+3+nBond*2-1)]
        numLine_off=[int(numLine_bond[1]+2),int(numLine_bond[1]+2+nOff*1-1)]
        numLine_ang=[int(numLine_off[1]+2 ),int(numLine_off[1]+2+nAng*1-1 )]
        numLine_tor=[int(numLine_ang[1]+2 ),int(numLine_ang[1]+2+nTor*1-1 )]
        numLine_hb =[int(numLine_tor[1]+2 ),int(numLine_tor[1]+2+nHB*1-1  )]
    with open(file_ffield,"r") as ffield:
        numLine=1
        gen_para=[]
        atom_type,atom_para=[],[]
        bond_type,bond_para=[],[]
        off_type,off_para=[],[]
        ang_type,ang_para=[],[]
        tor_type,tor_para=[],[]
        hb_type,hb_para=[],[]

        gen_lines = ffield.readlines()[numLine_gen[0]-1:numLine_gen[1]]
        for line in gen_lines:
            gen_para.append(line.split("!")[0])
        ffield.seek(0)
        atom_lines = ffield.readlines()[numLine_atom[0]-1:numLine_atom[1]]
        for line in atom_lines:
            if len(line.split())==9:
                atom_type.append(line.split()[0])
                atom_para.append(line.split()[1:])
            else:
                for para in line.split():
                    atom_para[len(atom_type)-1].append(para)
        ffield.seek(0)
        bond_lines = ffield.readlines()[numLine_bond[0]-1:numLine_bond[1]]
        for line in bond_lines:
            if len(line.split())==10:
                bond_type.append([int(line.split()[0]),int(line.split()[1])])
                bond_para.append(line.split()[2:])
            else:
                for para in line.split():
                    bond_para[len(bond_type)-1].append(para)
        ffield.seek(0)
        off_lines = ffield.readlines()[numLine_off[0]-1:numLine_off[1]]
        for line in off_lines:
            off_type.append([int(line.split()[0]),int(line.split()[1])])
            off_para.append(line.split()[2:])
        ffield.seek(0)
        ang_lines = ffield.readlines()[numLine_ang[0]-1:numLine_ang[1]]
        for line in ang_lines:
            ang_type.append([int(line.split()[0]),int(line.split()[1]),int(line.split()[2])])
            ang_para.append(line.split()[3:])
        ffield.seek(0)
        tor_lines = ffield.readlines()[numLine_tor[0]-1:numLine_tor[1]]
        for line in tor_lines:
            tor_type.append([int(line.split()[0]),int(line.split()[1]),int(line.split()[2]),int(line.split()[3])])
            tor_para.append(line.split()[4:])
        ffield.seek(0)
        hb_lines = ffield.readlines()[numLine_hb[0]-1:numLine_hb[1]]
        for line in hb_lines:
            hb_type.append([int(line.split()[0]),int(line.split()[1]),int(line.split()[2])])
            hb_para.append(line.split()[3:])
    general=[gen_type,gen_para]
    atom   =[atom_type,atom_para]
    bond   =[bond_type,bond_para]
    angle  =[off_type,off_para]
    off_diagonal=[ang_type,ang_para]
    torsion=[tor_type,tor_para]
    hbond  =[hb_type,hb_para]

    return general,atom,bond,angle,off_diagonal,torsion,hbond


def aName(int,pyffield):
    return pyffield[1][0][int-1]
def aNumber(str,pyffield):
    return pyffield[1][0].index(str)+1

def list_avail(pyffield,atom_type):
    dict_para={"general":0,"atom":1,"bond":2,"off":3,"angle":4,"torsion":5,"hbond":6}
    parameters=["general","atom","bond","off","angle","torsion","hbond"]
    ret_para_type=[[],[],[],[],[],[]]
    ret_para=[[],[],[],[],[],[]]
    atoms=[atom_type]
    numAtoms=[int(pyffield[1][0].index(atom))+1 for atom in atoms]
    for i in parameters:
        if dict_para[i]==1:
            count,j=0,0
            for atomic in pyffield[dict_para[i]][0]:
                if all(item in atomic for item in atoms):
                    ret_para_type[0].append(atomic)
                    ret_para[0].append(pyffield[dict_para[i]][1][j])
                    count+=1
                j+=1
        elif dict_para[i]==2:
            count,j=0,0
            for bond in pyffield[dict_para[i]][0]:
                if all(item in bond for item in numAtoms):
                    ret_para_type[1].append(aName(bond[0],pyffield)+"-"+aName(bond[1],pyffield))
                    ret_para[1].append(pyffield[dict_para[i]][1][j])
                    count+=1
                j+=1
        elif dict_para[i]==3:
            count,j=0,0
            for off in pyffield[dict_para[i]][0]:
                if all(item in off for item in numAtoms):
                    ret_para_type[2].append(aName(off[0],pyffield)+"-"+aName(off[1],pyffield))
                    ret_para[2].append(pyffield[dict_para[i]][1][j])
                    count+=1
                j+=1
        elif dict_para[i]==4:
            count,j=0,0
            for angle in pyffield[dict_para[i]][0]:
                if all(item in angle for item in numAtoms):
                    ret_para_type[3].append(aName(angle[0],pyffield)+"-"+aName(angle[1],pyffield)+"-"+aName(angle[2],pyffield))
                    ret_para[3].append(pyffield[dict_para[i]][1][j])
                    count+=1
                j+=1
        elif dict_para[i]==5:
            count,j=0,0
            for tor in pyffield[dict_para[i]][0]:
                if all(item in tor for item in numAtoms):
                    ret_para_type[4].append(aName(tor[0],pyffield)+"-"+aName(tor[1],pyffield)+"-"+aName(tor[2],pyffield)+"-"+aName(tor[3],pyffield))
                    ret_para[4].append(pyffield[dict_para[i]][1][j])
                    count+=1
                j+=1
        elif dict_para[i]==6:
            count,j=0,0
            for hb in pyffield[dict_para[i]][0]:
                if all(item in hb for item in numAtoms):
                    ret_para_type[5].append(aName(hb[0],pyffield)+"-"+aName(hb[1],pyffield)+"-"+aName(hb[2],pyffield))
                    ret_para[5].append(pyffield[dict_para[i]][1][j])
                    count+=1
                j+=1
    return ret_para_type, ret_para

def list_transfer(pyffield_get,pyffield_set,tr_atom):
    merge_avail=list_avail(pyffield_get,tr_atom)
    merge_atom=merge_avail[0][0][0]
    atoms_set=pyffield_set[1][0]
    filter_para_type=[[],[],[],[],[],[]]
    filter_para=[[],[],[],[],[],[]]
    filter_para_type[0].append(merge_avail[0][0][0])
    filter_para[0].append(merge_avail[1][0][0])
    count=0
    
    for bond in merge_avail[0][1]:
        para = [i for i in bond.split("-") if i != merge_atom]
        if all(item in atoms_set for item in para):
            filter_para_type[1].append(bond)
            filter_para[1].append(merge_avail[1][1][count])
        count+=1
    count=0
    for off in merge_avail[0][2]:
        para = [i for i in off.split("-") if i != merge_atom]
        if all(item in atoms_set for item in para):
            filter_para_type[2].append(off)
            filter_para[2].append(merge_avail[1][2][count])
        count+=1
    count=0
    for angle in merge_avail[0][3]:
        para = [i for i in angle.split("-") if i != merge_atom]
        if all(item in atoms_set for item in para):
            filter_para_type[3].append(angle)
            filter_para[3].append(merge_avail[1][3][count])
        count+=1
    count=0
    for tor in merge_avail[0][4]:
        para = [i for i in tor.split("-") if i != merge_atom]
        if all(item in atoms_set for item in para):
            filter_para_type[4].append(tor)
            filter_para[4].append(merge_avail[1][4][count])
        count+=1
    count=0
    for hb in merge_avail[0][5]:
        para = [i for i in hb.split("-") if i != merge_atom]
        if all(item in atoms_set for item in para):
            filter_para_type[5].append(hb)
            filter_para[5].append(merge_avail[1][5][count])
        count+=1
    return filter_para_type, filter_para

def format_merge(transferPara,pyffield_set):
    for i in range(6):
        if i ==0:
            txt_atom=""
            for j in range(4):
                if j==0:
                    txt_atom += " "+transferPara[0][i][0].ljust(2," ")
                else:
                    txt_atom += '   '
                for k in range(8*j,8*(j+1)):
                    txt_atom += '%9.4f' %float(transferPara[1][i][0][k])
                txt_atom+='\n'
        elif i ==1:
            txt_bond=""
            count=0
            for bond in transferPara[0][i]:            
                list_atoms_name=bond.split("-")
                list_atoms_number=[]
                for atom in list_atoms_name:
                    if atom==transferPara[0][0][0]:
                        list_atoms_number.append(len(pyffield_set[1][0])+1)
                    else:
                        list_atoms_number.append(aNumber(atom,pyffield_set))
                for j in range(2):
                    if j==0:
                        txt_bond += '%3d%3d' % (list_atoms_number[0],list_atoms_number[1])
                    else:
                        txt_bond +='      '
                    for k in range(8*j,8*(j+1)):
                        txt_bond += '%9.4f' %float(transferPara[1][i][count][k])
                    txt_bond+='\n'
                count+=1
        elif i ==2:
            txt_off=""
            count=0
            for off in transferPara[0][i]:   
                list_atoms_name=off.split("-")
                list_atoms_number=[]
                for atom in list_atoms_name:
                    if atom==transferPara[0][0][0]:
                        list_atoms_number.append(len(pyffield_set[1][0])+1)
                    else:
                        list_atoms_number.append(aNumber(atom,pyffield_set))
                txt_off += '%3d%3d' % (list_atoms_number[0],list_atoms_number[1])
                for k in range(0,6):
                    txt_off += '%9.4f' %float(transferPara[1][i][count][k])
                txt_off+='\n'
                count+=1        
        elif i ==3:
            txt_angle=""
            count=0
            for angle in transferPara[0][i]:   
                list_atoms_name=angle.split("-")
                list_atoms_number=[]
                for atom in list_atoms_name:
                    if atom==transferPara[0][0][0]:
                        list_atoms_number.append(len(pyffield_set[1][0])+1)
                    else:
                        list_atoms_number.append(aNumber(atom,pyffield_set))
                txt_angle += '%3d%3d%3d' % (list_atoms_number[0],list_atoms_number[1],list_atoms_number[2])
                for k in range(0,7):
                    txt_angle += '%9.4f' %float(transferPara[1][i][count][k])
                txt_angle+='\n'
                count+=1
        elif i ==4:
            txt_tor=""
            count=0
            for tor in transferPara[0][i]:
                list_atoms_name=tor.split("-")
                list_atoms_number=[]
                for atom in list_atoms_name:
                    if atom==transferPara[0][0][0]:
                        list_atoms_number.append(len(pyffield_set[1][0])+1)
                    else:
                        list_atoms_number.append(aNumber(atom,pyffield_set))
                txt_tor += '%3d%3d%3d%3d' % (list_atoms_number[0],list_atoms_number[1],list_atoms_number[2],list_atoms_number[3])
                for k in range(0,7):
                    txt_tor += '%9.4f' %float(transferPara[1][i][count][k])
                txt_tor+='\n'
                count+=1
        elif i ==5:
            txt_hb=""
            count=0
            for hb in transferPara[0][i]:
                list_atoms_name=hb.split("-")
                list_atoms_number=[]
                for atom in list_atoms_name:
                    if atom==transferPara[0][0][0]:
                        list_atoms_number.append(len(pyffield_set[1][0])+1)
                    else:
                        list_atoms_number.append(aNumber(atom,pyffield_set))
                txt_hb += '%3d%3d%3d' % (list_atoms_number[0],list_atoms_number[1],list_atoms_number[2])
                for k in range(0,4):
                    txt_hb += '%9.4f' %float(transferPara[1][i][count][k])
                txt_hb+='\n'
                count+=1
    return txt_atom,txt_bond,txt_off,txt_angle,txt_tor,txt_hb

def format_ffield(pyffield_para,p_type):
    if p_type=="gen":
        txt=""
        count=0
        for gen in pyffield_para[1]:
            txt+="%10.4f"%float(pyffield_para[1][count])+" !"+pyffield_para[0][count]+"\n"
            count+=1
    elif p_type=="atom":
        txt=""
        count=0
        for atom in pyffield_para[0]:
            for j in range(4):
                if j==0:

                    txt += ' '+pyffield_para[0][count].ljust(2," ")
                else:
                    txt += '   '
                for k in range(8*j,8*(j+1)):
                    txt += '%9.4f' %float(pyffield_para[1][count][k])
                txt+='\n'
            count+=1
    elif p_type=="bond":
        txt=""
        count=0
        for bond in pyffield_para[0]:            
            list_atoms_number=[int(i) for i in bond]
            for j in range(2):
                if j==0:
                    txt += '%3d%3d' % (list_atoms_number[0],list_atoms_number[1])
                else:
                    txt +='      '
                for k in range(8*j,8*(j+1)):
                    txt += '%9.4f' %float(pyffield_para[1][count][k])
                txt+='\n'
            count+=1

    elif p_type=="off":
        txt=""
        count=0
        for off in pyffield_para[0]:            
            list_atoms_number=[int(i) for i in off]
            txt += '%3d%3d' % (list_atoms_number[0],list_atoms_number[1])
            for k in range(0,6):
                txt += '%9.4f' %float(pyffield_para[1][count][k])
            txt+='\n'
            count+=1

    elif p_type=="angle":
        txt=""
        count=0
        for angle in pyffield_para[0]:            
            list_atoms_number=[int(i) for i in angle]
            txt += '%3d%3d%3d' % (list_atoms_number[0],list_atoms_number[1],list_atoms_number[2])
            for k in range(0,7):
                txt += '%9.4f' %float(pyffield_para[1][count][k])
            txt+='\n'
            count+=1

    elif p_type=="tor":
        txt=""
        count=0
        for tor in pyffield_para[0]:            
            list_atoms_number=[int(i) for i in tor]
            txt += '%3d%3d%3d%3d' % (list_atoms_number[0],list_atoms_number[1],list_atoms_number[2],\
                                     list_atoms_number[3])
            for k in range(0,7):
                txt += '%9.4f' %float(pyffield_para[1][count][k])
            txt+='\n'
            count+=1

    elif p_type=="hb":
        txt=""
        count=0
        for hb in pyffield_para[0]:            
            list_atoms_number=[int(i) for i in hb]
            txt += '%3d%3d%3d' % (list_atoms_number[0],list_atoms_number[1],list_atoms_number[2])
            for k in range(0,4):
                txt += '%9.4f' %float(pyffield_para[1][count][k])
            txt+='\n'
            count+=1
    return txt

def writeMerge(toAppend,pyffield_set,ffield_merge,atom):
    with open(ffield_merge,"w") as fout:
        para_atom,para_bond,para_off,para_angle,para_tor,para_hb=toAppend

        comment="This is the merged force field generated by taking W atom from a force field to another"
        com_gen ="%3d"%int(len(pyffield_set[0][0]))+"       ! Number of general parameters"
        com_atom="%3d"%int(len(pyffield_set[1][0])+(sum('\n' in s for s in para_atom))/4)+"    ! Nr of atoms; atomID;ro(sigma); Val;atom mass;Rvdw;Dij;gamma\n"\
                "            alfa;gamma(w);Val(angle);p(ovun5);n.u.;chiEEM;etaEEM;n.u.\n"\
                "            ro(pipi);p(lp2);Heat increment;p(boc4);p(boc3);p(boc5),n.u.;n.u.\n"\
                "            p(ovun2);p(val3);n.u.;Val(boc);p(val5);n.u.;n.u.;n.u."
        com_bond="%3d"%int(len(pyffield_set[2][0])+(sum('\n' in s for s in para_bond))/2)+"      ! Nr of bonds; at1;at2;De(sigma);De(pi);De(pipi);p(be1);p(b\n"\
                "                      p(be2);p(bo3);p(bo4);n.u.;p(bo1);p(bo2)"
        com_off ="%3d"%int(len(pyffield_set[3][0])+(sum('\n' in s for s in para_off))/1)+"    ! Nr of off-diagonal terms. at1;at2;Dij;RvdW;alfa;ro(sigma);r"
        com_angle ="%3d"%int(len(pyffield_set[4][0])+(sum('\n' in s for s in para_angle))/1)+"    ! Nr of angles. at1;at2;at3;Thetao,o;p(val1);p(val2);p(coa1);"
        com_tor ="%3d"%int(len(pyffield_set[5][0])+(sum('\n' in s for s in para_tor))/1)+"    ! Nr of torsions. at1;at2;at3;at4;;V1;V2;V3;p(tor1);p(cot1);n"
        com_hb ="%3d"%int(len(pyffield_set[6][0])+(sum('\n' in s for s in para_hb))/1)+"    ! Nr of hydrogen bonds. at1;at2;at3;r(hb);p(hb1);p(hb2);p(hb3"
        fout.write(comment+"\n")
        fout.write(com_gen+"\n")
        fout.write(format_ffield(pyffield_set[0],"gen"))
        fout.write(com_atom+"\n")
        fout.write(format_ffield(pyffield_set[1],"atom")+para_atom)
        fout.write(com_bond+"\n")
        fout.write(format_ffield(pyffield_set[2],"bond")+para_bond)
        fout.write(com_off+"\n")
        fout.write(format_ffield(pyffield_set[3],"off") +para_off )
        fout.write(com_angle+"\n")
        fout.write(format_ffield(pyffield_set[4],"angle")+para_angle)
        fout.write(com_tor+"\n")
        fout.write(format_ffield(pyffield_set[5],"tor") +para_tor )
        fout.write(com_hb+"\n")
        fout.write(format_ffield(pyffield_set[6],"hb")  +para_hb)

#file="ffield_get"
filename=input("Enter the source force field: ")
#print(file)
#atom="W"
atom=input("Which atom to merge: ")
pyffield_get=read_ffield(filename)
para_avail=list_avail(pyffield_get,atom)
dict_para={"general":0,"atom":1,"bond":2,"off":3,"angle":4,"torsion":5,"hbond":6}
parameters=["General","Atom","Bond","Off-diagonal","Angle","Torsion","H-bond"]    

file_set=input("Enter the destination force field: ")
#file_set="ffield_set"
pyffield_set=read_ffield(file_set)
merge_avail=list_transfer(pyffield_get,pyffield_set,atom)

print("====================== %I am merging these parameters% ======================")


for i in range(1,7):
    para_list=""
    print("\n")
    print(parameters[i]+": \n|")
    if len(merge_avail[0][i-1])==0:
        print("Not Available")
    else:
        for para in merge_avail[0][i-1]:
            para_list+=para+" | "
    print(para_list)
print("============================================================================\n")
#proceed=input("Do you want to continue? y/n ")
proceed="y"
if proceed=="n":
    pass
elif proceed=="y":
    toAppend=format_merge(merge_avail,pyffield_set)
    print("\n============== %Here is how the formatted data looks like% ==============")
    for i in range(0,5):
        print(parameters[i+1]+": \n|")
        print(toAppend[i])
    #nameMerge=input("Enter the name of the merged force field: ")
    nameMerge="ffield_merged"

    writeMerge(toAppend,pyffield_set,nameMerge,atom)
    print("Merged SuccessFully!\nNew force field is saved as ffield_merged.\nYou're Welcome!")
