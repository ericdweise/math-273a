import matplotlib.pyplot as plt

def euler_init(h, n_points):
    grid = [ h*v for v in range(-n_points,n_points) ]
    graph = [ 1 if x > 0 else 0 for x in grid ]
    return grid, graph

def euler_left(k, h, prev_graph):
    new_graph = [-1]*len(prev_graph) 
    for i in range(1, len(prev_graph)):
        new_graph[i] = prev_graph[i] - k*(prev_graph[i] - prev_graph[i-1])/h

    # assume left most point has same value its left.
    new_graph[0] = prev_graph[0]
    return new_graph

def euler_right(k, h, prev_graph):
    new_graph = [-1]*len(prev_graph) 
    for i in range(0, len(prev_graph)-1):
        new_graph[i] = prev_graph[i] - k*(prev_graph[i+1] - prev_graph[i])/h

    # assume right most point has same value to its right.
    new_graph[-1] = prev_graph[-1]
    return new_graph

def euler(k, h, prev_graph):
    new_graph = [-1]*len(prev_graph) 
    for i in range(1, len(prev_graph)-1):
        new_graph[i] = prev_graph[i] - k*(prev_graph[i+1] - prev_graph[i-1])/2/h

    # assume left most point has same value its left.
    new_graph[0] = prev_graph[0] - k*(prev_graph[1] - prev_graph[0])/2/h

    # assume right most point has same value to its right.
    new_graph[-1] = prev_graph[-1] - k*(prev_graph[-1] - prev_graph[-2])/2/h

    return new_graph


# Initialize values
h = 0.05
k = 0.05
n = 100
grid, original_graph = euler_init(h,n)
plt.plot(grid, original_graph)
plt.savefig('original_graph.png')
plt.clf()


##########
# PART B #
##########
prev_graph = original_graph
for n_iter in range(10):
    next_graph = euler_left(k, h, prev_graph)
    plt.plot(grid, next_graph)
    plt.savefig('left_{}.png'.format(n_iter+1))
    plt.clf()
    prev_graph = next_graph


##########
# PART C #
##########
prev_graph = original_graph
for n_iter in range(10):
    next_graph = euler_right(k, h, prev_graph)
    plt.plot(grid, next_graph)
    plt.savefig('right_{}.png'.format(n_iter+1))
    plt.clf()
    prev_graph = next_graph


##########
# PART D #
##########
prev_graph = original_graph
for n_iter in range(10):
    next_graph = euler(k, h, prev_graph)
    plt.plot(grid, next_graph)
    plt.savefig('euler_{}.png'.format(n_iter+1))
    plt.clf()
    prev_graph = next_graph
