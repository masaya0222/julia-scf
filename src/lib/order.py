import constant
from constant import hash, io, MAX_AM
from mem_man import get_total_memory, get_mem, init_mem, free_mem
import numpy  as np

MAX_NODE = 20000
NONODE = -1000000

last_hrr_node = 0
last_vrr_node = 0
first_hrr_to_compute = 0
first_vrr_to_compute = 0

hrr_hash_table = np.zeros((MAX_AM+1,MAX_AM+1,MAX_AM+1,MAX_AM+1),dtype=int)
vrr_hash_table = np.zeros((MAX_AM+1,MAX_AM+1,2*MAX_AM+1),dtype=int)

class hrr_class:
    def __init__(self):
        self.A,self.B,self.C,self.D = (0,0,0,0)
        self.size = 0
        self.pointer = 0
        self.children = [0,0]
        self.parents_counter = 0
        self.num_parents = 0
        self.parents = [0,0,0,0,0]
        self.llink = 0
        self.rlink = 0
        self.marked = 0
        self.target = 0

class vrr_class:
    def __init__(self):
        self.A, self.C = (0,0)
        self.m = 0
        self.size = 0
        self.pointer = 0
        self.children = [0 for i in range(5)]
        self.parents_counter = 0
        self.num_parents = 0
        self.parents = [0 for i in range(9)]
        self.llink = 0
        self.rlink = 0
        self.marked = 0
        self.target = 0


def emit_order(libint_stack_size):
    new_am = constant.new_am
    opt_am = constant.opt_am
    am_letter = constant.am_letter
    lib = constant.lib
    global last_hrr_node, first_hrr_to_compute, last_vrr_node, first_vrr_to_compute
    global hrr_hash_table, vrr_hash_table
    hrr_nodes = [hrr_class() for i in range(MAX_NODE)]
    vrr_nodes = [vrr_class() for i in range(MAX_NODE)]
    target_vrr_nodes = np.zeros(MAX_AM*MAX_AM,dtype=int)

    max_stack_size = 0
    current_highest_am = 0
    hrr_order = "hrr_order.jl"
    hrr_order_code = open(lib+"/"+hrr_order, 'w')

    for la in range(new_am+1):
        lb_max = la//2
        lb_min = 0 if la <= new_am//2 else la - new_am//2
        lc_min = 1 if la == 0 else la
        for lb in range(lb_min,lb_max+1):
            for lc in range(lc_min, new_am+1):
                ld_max = lc//2
                ld_min = 0 if lc <= new_am//2 else lc - new_am//2
                for ld in range(ld_min, ld_max+1):
                    current_highest_am = max(la-lb, lb, lc-ld, ld)
                    hrr_function_name = "hrr_order_"+am_letter[la-lb]+am_letter[lb]+am_letter[lc-ld]+am_letter[ld]
                    vrr_function_name = "vrr_order_"+am_letter[la-lb]+am_letter[lb]+am_letter[lc-ld]+am_letter[ld]

                    hrr_code_name = hrr_function_name+".jl"
                    vrr_code_name = vrr_function_name+".jl"
                    hrr_function_name += "!"                    
                    hrr_code = open(lib+"/" + hrr_code_name, 'w')
                    vrr_code = open(lib+"/"+vrr_code_name, 'w')
                    
                    hrr_order_code.write('include("%s")\n'%hrr_code_name)

                    hrr_code.write('include("libints_header.jl")\n')
                    hrr_code.write('include("hrr_header.jl")\n')
                    hrr_code.write('include("%s")\n\n'%vrr_code_name)
                    hrr_code.write("# Computes quartets of (%s%s|%s%s) integrals\n\n"%(am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]))

                    hrr_code.write("function %s(Libint::Libint_t, num_prim_comb::Int)\n"%(hrr_function_name))
                    hrr_code.write("\tData = Libint.PrimQuartet\n")
                    hrr_code.write("\tData_ind = 1\n")
                    hrr_code.write("\tint_stack = Libint.int_stack\n")
                    #hrr_code.write("\tint_stack_ind = 1\n\n")
                    
                    vrr_code.write('include("libints_header.jl")\n')
                    vrr_code.write('include("vrr_header.jl")\n\n')
                    vrr_code.write("# Computes quartets necessary to compute (%s%s|%s%s) integrals\n\n"%(am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld]))
                    vrr_code.write("function %s!(Libint::Libint_t, Data::prim_data)\n"%(vrr_function_name))
                    
                    hrr_hash_table = np.zeros((MAX_AM+1,MAX_AM+1,MAX_AM+1,MAX_AM+1),dtype=int)
                    vrr_hash_table = np.zeros((MAX_AM+1,MAX_AM+1,2*MAX_AM+1),dtype=int)

                    hrr_nodes[0].A = la-lb
                    hrr_nodes[0].B = lb
                    hrr_nodes[0].C = lc-ld
                    hrr_nodes[0].D = ld
                    hrr_nodes[0].llink = -1
                    hrr_nodes[0].rlink = -1
                    first_hrr_to_compute = 0
                    last_hrr_node = 0
                    mk_hrr_node(hrr_nodes[0], hrr_nodes, 0)
                    hrr_nodes[0].target = 1

                    for k in range(2):
                        if (hrr_nodes[0].children[k] > 0):
                            mark_hrr_parents(hrr_nodes[0].children[k], hrr_nodes, 0)
                    
                    init_mem(1)
                    
                    for i in range(last_hrr_node-1, -1, -1):
                        if hrr_nodes[i].B == 0 and hrr_nodes[i].D == 0:
                            hrr_nodes[i].marked = 1
                            hrr_nodes[i].pointer = get_mem(hrr_nodes[i].size)
                            hrr_code.write("\tLibint.vrr_classes[%d,%d] = %d\n"%(hrr_nodes[i].A+1, hrr_nodes[i].C+1, hrr_nodes[i].pointer))

                    base_mem = get_total_memory()
                    hrr_code.write("\tfor i = 1:%d\n"%base_mem)
                    hrr_code.write("\t\tint_stack[i] = 0\n")
                    hrr_code.write("\tend\n\n")
                    hrr_code.write("\tLibint.vrr_stack = %d\n"%base_mem)
                    
                    target_data = alloc_mem_hrr(hrr_nodes)
                    hrr_mem = get_total_memory()
                    if (max_stack_size < hrr_mem):
                        max_stack_size = hrr_mem
                    hrr_code.write("\tfor i = 1:num_prim_comb\n")
                    hrr_code.write("\t\tvrr_order_%s%s%s%s!(Libint, Data[Data_ind])\n"%((am_letter[la-lb],am_letter[lb],am_letter[lc-ld],am_letter[ld])))
                    hrr_code.write("\t\tData_ind += 1\n")
                    hrr_code.write("\tend\n")
                    j = first_hrr_to_compute
                    
                    while True:
                        hrr_code.write("# --- compute (%s%s|%s%s) ---\n"%(am_letter[hrr_nodes[j].A],am_letter[hrr_nodes[j].B], am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D]))
                        if (hrr_nodes[j].B == 0 and hrr_nodes[j].D != 0):
                            hrr_code.write("\thrr3_build_%s%s!(Libint.CD, int_stack, %d, int_stack, %d,"%(am_letter[hrr_nodes[j].C], am_letter[hrr_nodes[j].D],\
                                hrr_nodes[j].pointer, hrr_nodes[hrr_nodes[j].children[0]].pointer))
                            hrr_code.write(" int_stack, %d, %d)\n"%(hrr_nodes[hrr_nodes[j].children[1]].pointer, io(1+hrr_nodes[j].A)*io(1+hrr_nodes[j].B)))
                        elif hrr_nodes[j].B != 0:
                            hrr_code.write("\thrr1_build_%s%s!(Libint.AB, int_stack, %d, int_stack, %d,"%(am_letter[hrr_nodes[j].A], am_letter[hrr_nodes[j].B],\
                                hrr_nodes[j].pointer, hrr_nodes[hrr_nodes[j].children[0]].pointer))
                            hrr_code.write(" int_stack, %d, %d)\n"%(hrr_nodes[hrr_nodes[j].children[1]].pointer, io(1+hrr_nodes[j].C)*io(1+hrr_nodes[j].D)))
                        j = hrr_nodes[j].rlink
                        if j == -1:
                            break
                    hrr_code.write("\treturn %d\nend\n"%target_data)
                    hrr_code.close()
                    print("Done with %s"%hrr_code_name)
                    for i in range(last_hrr_node):
                        hrr_nodes[i].llink = 0
                        hrr_nodes[i].rlink = 0
                    
                    last_vrr_node = 0
                    num_vrr_targets = 0
                    for i in range(last_hrr_node):
                        if (hrr_nodes[i].B == 0 and hrr_nodes[i].D == 0):
                            j = vrr_hash_table[hrr_nodes[i].A][hrr_nodes[i].C][0]
                            if (j == 0):
                                target_vrr_nodes[num_vrr_targets] = last_vrr_node
                                vrr_nodes[last_vrr_node].A = hrr_nodes[i].A
                                vrr_nodes[last_vrr_node].C = hrr_nodes[i].C
                                vrr_nodes[last_vrr_node].m = 0
                                vrr_nodes[last_vrr_node].llink = -1
                                vrr_nodes[last_vrr_node].rlink = -1
                                if num_vrr_targets:
                                    vrr_nodes[target_vrr_nodes[num_vrr_targets-1]].llink = last_vrr_node
                                    vrr_nodes[last_vrr_node].rlink = target_vrr_nodes[num_vrr_targets-1]
                                    vrr_nodes[last_vrr_node].llink = -1
                                else:
                                    vrr_nodes[last_vrr_node].rlink = -1
                                    vrr_nodes[last_vrr_node].llink = -1
                                first_vrr_to_compute = last_vrr_node
                                mk_vrr_node(vrr_nodes[last_vrr_node], vrr_nodes, 0)
                                vrr_nodes[first_vrr_to_compute].target = 1
                                num_vrr_targets+= 1
                            else:
                                vrr_nodes[j-1].target = 1
                            if (first_vrr_to_compute == last_vrr_node and i == last_hrr_node-1):
                                print("Edward, you fucked up")
                    
                    for i in range(num_vrr_targets):
                        j = target_vrr_nodes[i]
                        for k in range(5):
                            if vrr_nodes[j].children[k] > 0:
                                mark_vrr_parents(vrr_nodes[j].children[k], vrr_nodes, j)
                    init_mem(1)

                    target_data = alloc_mem_vrr(vrr_nodes)
                    vrr_mem = base_mem + get_total_memory()
                    if (max_stack_size < vrr_mem):
                        max_stack_size = vrr_mem
                    vrr_code.write("\tvrr_stack = Libint.vrr_stack\n")

                    j = first_vrr_to_compute
                    while True:
                        if vrr_nodes[j].A <= opt_am and vrr_nodes[j].C <= opt_am:
                            vrr_code.write("\tbuild_%s0%s0!(Data, Libint.int_stack, "%(am_letter[vrr_nodes[j].A],am_letter[vrr_nodes[j].C]))
                        else:
                            vrr_code.write(" am[0] = %d; am[1] = %d\n"%(vrr_nodes[j].A, vrr_nodes[j].C))
                            vrr_code.write("vrr_build_xxxx(am,Data,")
                        vrr_code.write("vrr_stack+%d"%vrr_nodes[j].pointer)
                        for k in range(5):
                            if (vrr_nodes[j].children[k] > 0):
                                vrr_code.write(", Libint.int_stack, vrr_stack+%d"%vrr_nodes[vrr_nodes[j].children[k]].pointer)
                            elif vrr_nodes[j].children[k] == NONODE:
                                vrr_code.write(", [0.0], -1")
                            else:
                                vrr_code.write(", Data.F, %d"%((-1)*vrr_nodes[j].children[k]))
                        vrr_code.write(")\n")
                        if (vrr_nodes[j].target == 1):
                            vrr_code.write("\ttmp = vrr_stack + %d\n"%vrr_nodes[j].pointer)
                            vrr_code.write("\ttarget_ptr = Libint.vrr_classes[%d,%d]\n"%(vrr_nodes[j].A+1, vrr_nodes[j].C+1))
                            vrr_code.write("\tfor i = 1:%d\n"%((vrr_nodes[j].A+1)*(vrr_nodes[j].A+2)*(vrr_nodes[j].C+1)*(vrr_nodes[j].C+2)//4))
                            vrr_code.write("\t\tLibint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]\n")
                            vrr_code.write("\tend\n")
                        j = vrr_nodes[j].rlink
                        if j == -1:
                            break
                    vrr_code.write("end\n\n")
                    
                    vrr_code.close()
                    #print("Done with %s"%vrr_code_name)
                    for i in range(last_vrr_node):
                        vrr_nodes[i].llink = 0
                        vrr_nodes[j].rlink = 0

                    if libint_stack_size[current_highest_am] < max_stack_size:
                        libint_stack_size[current_highest_am] = max_stack_size
                    
                    max_stack_size = 0

    hrr_order_code.close()

# Recursive function that build the HRR subgraph given the parent
def mk_hrr_node(node, allnodes, new):
    global last_hrr_node
    global hrr_hash_table
    O = [hrr_class(), hrr_class()]
    made = 0
    # Search for the parent node on stack
    #If it's not there - we'll add it to the end of the stack 
    thisnode = last_hrr_node
    # it's already placed on the stack allnodes - make sure children don't get created again (made = 1) 
    if (hrr_hash_table[node.A][node.B][node.C][node.D]):
        i = hrr_hash_table[node.A][node.B][node.C][node.D] - 1
        thisnode = i
        made = 1
    
    # it's not computed, add it, and make it the first to compute! 
    if not made:
        allnodes[thisnode].A = node.A
        allnodes[thisnode].B = node.B
        allnodes[thisnode].C = node.C
        allnodes[thisnode].D = node.D
        hrr_hash_table[node.A][node.B][node.C][node.D] = thisnode + 1
        allnodes[thisnode].num_parents = 0
        allnodes[thisnode].parents_counter = 0
        allnodes[thisnode].marked = 0
        allnodes[thisnode].pointer = 0
        for ind in range(5):
            allnodes[thisnode].parents[ind] = 0
        allnodes[thisnode].children[0] = NONODE
        allnodes[thisnode].children[1] = NONODE
        allnodes[thisnode].size = io(1+node.A)*io(1+node.B)*io(1+node.C)*io(1+node.D)
        allnodes[thisnode].target = 0
        # We just added a node ..
        last_hrr_node += 1
        if (last_hrr_node == MAX_NODE):
            print("Maximum stack size is reached. Change MAXNODE and recompile")
            exit(1)
    # If the parent class wasn't on stack already (!new) - increase the parent counter
    if(not new):
        allnodes[thisnode].num_parents += 1
        allnodes[thisnode].parents_counter += 1
    
    # now make all child nodes
    if (not made):
        if (node.B):
            O[0].A = node.A+1
            O[0].B = node.B-1
            O[0].C = node.C
            O[0].D = node.D
            allnodes[thisnode].children[0] = mk_hrr_node(O[0], allnodes, made)
            O[1].A = node.A
            O[1].B = node.B-1
            O[1].C = node.C
            O[1].D = node.D
            allnodes[thisnode].children[1] = mk_hrr_node(O[1], allnodes, made)
        elif (node.D):
            O[0].A = node.A
            O[0].B = node.B
            O[0].C = node.C+1
            O[0].D = node.D-1
            allnodes[thisnode].children[0] = mk_hrr_node(O[0], allnodes, made)
            O[1].A = node.A
            O[1].B = node.B
            O[1].C = node.C
            O[1].D = node.D-1
            allnodes[thisnode].children[1] = mk_hrr_node(O[1], allnodes, made)
    return thisnode

# Recursive function that build the VRR subgraph given the parent
def mk_vrr_node(node, allnodes, new):
    global last_vrr_node, vrr_hash_table
    O = [vrr_class() for i in range(5)]
    made = 0
    if node.A + node.C ==0:
        return -1*node.m
    thisnode = last_vrr_node
    if (vrr_hash_table[node.A][node.C][node.m]):
        i = vrr_hash_table[node.A][node.C][node.m] -1 
        thisnode = i
        made = 1
    
    if not made:
        allnodes[thisnode].A = node.A
        allnodes[thisnode].C = node.C
        allnodes[thisnode].m = node.m
        vrr_hash_table[node.A][node.C][node.m] = thisnode + 1
        allnodes[thisnode].num_parents = 0
        allnodes[thisnode].parents_counter = 0
        allnodes[thisnode].marked = 0
        allnodes[thisnode].target = 0
        allnodes[thisnode].pointer = 0
        for ind in range(9):
            allnodes[thisnode].parents[ind] = 0
        allnodes[thisnode].children[0] = NONODE
        allnodes[thisnode].children[1] = NONODE
        allnodes[thisnode].children[2] = NONODE
        allnodes[thisnode].children[3] = NONODE
        allnodes[thisnode].children[4] = NONODE
        allnodes[thisnode].size = io(1+node.A)*io(1+node.C)
        # We just added a node ..
        last_vrr_node += 1
        # If stack is overfull - exit 
        if(last_vrr_node==MAX_NODE):
            print(" Maximum stack size is reached. Change MAXNODE and recompile.")
            exit(1)
    
    # If the parent class wasn't on stack already (!new) - increase the parent counter 
    if(not new):
        allnodes[thisnode].num_parents += 1
        allnodes[thisnode].parents_counter += 1
    
    # now make all child nodes 
    if (not made):
        if(node.A):
            O[0].A = node.A-1
            O[0].C = node.C
            O[0].m = node.m
            allnodes[thisnode].children[0] = mk_vrr_node(O[0], allnodes, made)
            O[1].A = node.A-1
            O[1].C = node.C
            O[1].m = node.m+1
            allnodes[thisnode].children[1] = mk_vrr_node(O[1], allnodes, made)
            if(node.A>1):
                O[2].A = node.A-2
                O[2].C = node.C
                O[2].m = node.m
                allnodes[thisnode].children[2] = mk_vrr_node(O[2], allnodes, made)
                O[3].A = node.A-2
                O[3].C = node.C
                O[3].m = node.m+1
                allnodes[thisnode].children[3] = mk_vrr_node(O[3], allnodes, made)
            
            if(node.C):
                O[4].A = node.A-1
                O[4].C = node.C-1
                O[4].m = node.m+1
                allnodes[thisnode].children[4] = mk_vrr_node(O[4], allnodes, made)

        elif(node.C):
            O[0].A = node.A
            O[0].C = node.C-1
            O[0].m = node.m
            allnodes[thisnode].children[0] = mk_vrr_node(O[0], allnodes, made)
            O[1].A = node.A
            O[1].C = node.C-1
            O[1].m = node.m+1
            allnodes[thisnode].children[1] = mk_vrr_node(O[1], allnodes, made)
            if(node.C>1):
                O[2].A = node.A
                O[2].C = node.C-2
                O[2].m = node.m
                allnodes[thisnode].children[2] = mk_vrr_node(O[2], allnodes, made)
                O[3].A = node.A
                O[3].C = node.C-2
                O[3].m = node.m+1
                allnodes[thisnode].children[3] = mk_vrr_node(O[3], allnodes, made)
    return thisnode

# Make hrr_nodes[rent] a parent of hrr_nodes[n] and proceed recursively
def mark_hrr_parents(n, allnodes, rent):
    global first_hrr_to_compute
    # handle case where it's in the parent list already
    for i in range(allnodes[n].num_parents-1, allnodes[n].parents_counter-1,-1):
        if rent == allnodes[n].parents[i]:
            return
    # if the parent rent is not in the list - add it to the list!
    allnodes[n].parents_counter -= 1
    i = allnodes[n].parents_counter
    allnodes[n].parents[i] = rent
    # hits from all of the parents has been received - schedule it for computation and mark all of its children
    if i == 0 and (allnodes[n].B != 0 or allnodes[n].D != 0):
        allnodes[n].llink = -1
        allnodes[n].rlink = first_hrr_to_compute
        allnodes[first_hrr_to_compute].llink = n
        first_hrr_to_compute = n

        for i in range(2):
            if (allnodes[n].children[i]>0):
                mark_hrr_parents(allnodes[n].children[i], allnodes, n)
    return

# Make vrr_nodes[rent] a parent of vrr_nodes[n] and proceed recursively 
def mark_vrr_parents(n, allnodes, rent):
    global first_vrr_to_compute
    for i in range(allnodes[n].num_parents-1, allnodes[n].parents_counter-1,-1):
        if rent == allnodes[n].parents[i]:
            return
    allnodes[n].parents_counter -= 1
    i = allnodes[n].parents_counter
    allnodes[n].parents[i] = rent
    if i == 0:
        allnodes[n].llink = -1
        allnodes[n].rlink = first_vrr_to_compute
        allnodes[first_vrr_to_compute].llink = n
        first_vrr_to_compute = n
        for i in range(5):
            if (allnodes[n].children[i] > 0):
                mark_vrr_parents(allnodes[n].children[i], allnodes, n)
    return

def alloc_mem_hrr(nodes):
    global first_hrr_to_compute
    
    j = first_hrr_to_compute
    while True:
        if nodes[j].marked == 0:
            nodes[j].marked = 1
            nodes[j].pointer = get_mem(nodes[j].size)
        for k in range(2):
            child = nodes[j].children[k]
            if child > 0:
                if (nodes[child].target == 0):
                    free_it = 1
                    for l in range(nodes[child].num_parents):
                        if not nodes[nodes[child].parents[l]].marked:
                            free_it = 0
                    if free_it:
                        free_mem(nodes[child].pointer, nodes[child].size)
        j = nodes[j].rlink
        
        if j == -1:
            break
    return nodes[0].pointer

# This functions controls memory placement of computed classes on the CINTS stack 
def alloc_mem_vrr(nodes):
    global last_vrr_node, first_vrr_to_compute
    for i in range(last_vrr_node):
        nodes[i].marked = 0
    j = first_vrr_to_compute
    while True:
        nodes[j].marked = 1
        nodes[j].pointer = get_mem(nodes[j].size)

        for k in range(5):
            child = nodes[j].children[k]
            if child > 0:
                free_it = 1
                for l in range(nodes[child].num_parents):
                    if not nodes[nodes[child].parents[l]].marked:
                        free_it = 0
                if (free_it):
                    free_mem(nodes[child].pointer, nodes[child].size)
        j = nodes[j].rlink
        if j == -1:
            break
    return nodes[0].pointer

