def get_rand_BEC(struc, eqop, symops, max_charge=1):
    symstruc = sga(struc)
    maptoprim = sgastruc.get_symmetry_dataset()['mapping_to_primitive']
    eqatoms = sgastruc.get_symmetry_dataset()['equivalent_atoms']
    
    l = len(struc)
    BEC = np.zeros((l,3,3))
    primsites = []
    for i in range(l):
        good = 0    
        for j in reversed(range(i)):
            if maptoprim[j] == maptoprim[i]:
                BEC[i] = BEC[j]
                good = 1
                break
            if eqatoms[j] == eqatoms[i] and maptoprim[j] != maptoprim[i]:
                for k in range(len(eqop[i])):
                    if eqop[i][k][0] == j:
                        tempbec = BEC[j]
                        for m in range(eqop[i][k][2]):
                            tempbec = eqop[i][k][1].transform_tensor(tempbec)
                            good = 1

                        
                        BEC[i] = tempbec
                    if good == 1:
                        break
            if good == 1:
                break
        if good == 1:
            continue
                
                
        symsite = 0 
        temp_tensor = Tensor(np.random.rand(3,3)-0.5)         
        temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in symops[i]]) \
            /len(symops[i])
            
        BEC[i] = np.array(temp_tensor)

    count = 0

    
#     Enforce Acoustic Sum
    disp_charge = np.einsum("ijk->jk",BEC)/l
    
    add = np.zeros([l,3,3])
            
    for i in range(l):
        add[i] = disp_charge
        
        eq = 0
        for j in range(i):
            if maptoprim[i] == maptoprim[j] and i != j:
                add[i] = add[j]
                eq = 1
                break

        if eq == 1:
            continue

        good = 0    
        for j in reversed(range(i)):

            if eqatoms[j] == eqatoms[i]:
                for k in range(len(eqop[i])):
                    if eqop[i][k][0] == j:
                        tempbec = add[j]
                        for m in range(eqop[i][k][2]):
                            tempbec = eqop[i][k][1].transform_tensor(tempbec)
                            good = 1
                            break
                        add[i] = tempbec
                    if good == 1:
                        break
            if good == 1:
                break
        if good == 1:
            continue

        add[i] = sum([Tensor(add[i]).transform(symm_op) for symm_op in symops[i]]) / len(symops[i])
    
    BEC = BEC - add
    
    disp_charge = np.einsum("ijk->jk",BEC)/l
 
    return BEC*max_charge

def get_rand_IST(struc, eqop, symops, max_force=1):
    symstruc = sga(struc)
    maptoprim = sgastruc.get_symmetry_dataset()['mapping_to_primitive']
    eqatoms = sgastruc.get_symmetry_dataset()['equivalent_atoms']
    
    l = len(struc)
    IST = np.zeros((l,3,3,3))
    for i in range(l):
        good = 0    
        for j in reversed(range(i)):
            if maptoprim[j] == maptoprim[i]:
                IST[i] = IST[j]
                good = 1
                break
            if eqatoms[j] == eqatoms[i]:
                
                for k in range(len(eqop[i])):
                    if eqop[i][k][0] == j:
                        tempist = IST[j]
                        for m in range(eqop[i][k][2]):
                            tempist = eqop[i][k][1].transform_tensor(tempist)
                            good = 1
                        IST[i] = tempist
                    if good == 1:
                        break
            if good == 1:
                break
        if good == 1:
            continue
                
        temp_tensor = Tensor(np.random.rand(3,3,3)-0.5)    
        for m in range(3):
            for n in range(3):
                for o in range(n,3):
                    temp_tensor[m][o][n] = temp_tensor[m][n][o]
        temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in symops[i]]) \
            /len(symops[i])
 
        IST[i] = np.array(temp_tensor)
    count = 0

    return IST*max_force


def get_FCM_symmetry(fcm, operations, sharedops, max_force = 1):
    # set max force in reciprocal space
    numsites = int(len(fcm/3))
    D = fcm
    passed = []
    #symmetry operations
    for i in range(len(operations)):
        if type(operations[i][4]) == pymatgen.core.operations.SymmOp:
#             print(4, operations[i][0], operations[i][1])

            tempfcm = D[3*operations[i][2]:3*operations[i][2]+3, 3*operations[i][3]:3*operations[i][3]+3]
            tempfcm = operations[i][4].transform_tensor(tempfcm)

            D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = tempfcm
            D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = tempfcm.T

            continue

#         print(5, operations[i][0], operations[i][1])
        temp_tensor = Tensor(D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3])
        temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in \
                               sharedops[operations[i][0]][operations[i][1]]])
        temp_tensor_sum = temp_tensor_sum/(len(sharedops[operations[i][0]][operations[i][1]]))
        
        for pair in range(len(operations[i][5])):
            
            temp_tensor2 = temp_tensor_sum.T
            temp_tensor2 = operations[i][5][pair].transform_tensor(temp_tensor2)
            temp_tensor_sum = (temp_tensor_sum + temp_tensor2)/2

        D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = temp_tensor_sum
        D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = temp_tensor_sum.T

    return(D)