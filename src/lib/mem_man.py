import numpy as np
MAXALLOC = 100000

free_block = np.zeros(MAXALLOC, int)
block_length = np.zeros(MAXALLOC, int)
last_free = 0
max_mem = 0
mem_top = 0


def init_mem(memory):
    global last_free, free_block, block_length, max_mem, mem_top
    last_free = 1
    free_block = np.zeros(MAXALLOC, int)
    block_length = np.zeros(MAXALLOC, int)
    block_length[0] = memory
    max_mem = 0
    mem_top = memory

def add_mem(memory):
    global last_free, free_block, block_length, max_mem, mem_top
    addto = 0
    himem = 0
    mem_top += memory
    for i in range(last_free):
        if (free_block[i] > himem):
            himem = free_block[i]
            addto = i
    block_length[addto] = block_length[addto]+memory

def get_total_memory():
    global max_mem
    return max_mem

def get_mem(size):
    global last_free, free_block, block_length, max_mem, mem_top
    for i in range(last_free-1, -1, -1):
        if (block_length[i] == size):
            j = free_block[i]
            if (free_block[i]+block_length[i] == mem_top):
                add_mem(500)
            use(i, size)
            return j
    for i in range(last_free-1, -1, -1):
        if (block_length[i] > size):
            j = free_block[i]
            use(i, size)
            return j
    add_mem(1000)
    return get_mem(size)

def use(n, s):
    global last_free, free_block, block_length, max_mem, mem_top
    if (s == block_length[n]):
        for i in range(n,last_free):
            free_block[i] = free_block[i+1]
            block_length[i] = block_length[i+1]
        last_free -= 1
    elif (s<block_length[n]):
        free_block[n] = free_block[n] + s
        block_length[n] = block_length[n] -s
        if (free_block[n] > max_mem):
            max_mem = free_block[n]
    else:
        exit(1)

def free_mem(n, size):
    global last_free, free_block, block_length, max_mem, mem_top
    free_block[last_free] = n
    block_length[last_free] = size
    last_free += 1
    consolidate()

def consolidate():
    global last_free, free_block, block_length, max_mem, mem_top
    done = 1
    while True:
        done = 1
        for i in range(last_free):
            if last_free <= i:
                break
            right_bound_i = free_block[i]+block_length[i]
            for j in range(last_free):
                if last_free <= j:
                    break
                if free_block[j] == right_bound_i:
                    block_length[i] += block_length[j]
                    use(j, block_length[j])
                    done = 0
        if done:
            break
