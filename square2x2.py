import numpy


import pygalmesh


def square2x2():
    points = numpy.array([[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])
    constraints = [[0, 1], [1, 2], [2, 3], [3, 0]]
    mesh = pygalmesh.generate_2d(
        points,
        constraints,
        max_edge_size=0.3,
        num_lloyd_steps=10,
    )

    print(len(mesh.points), len( mesh.get_cells_type("triangle")))
    write_node(mesh.points)
    write_ele(len(mesh.points), mesh.get_cells_type("triangle"))
    mesh.write("square.svg")

def write_node(vertices):
    largo =len(vertices)
    f = open('input/square2x2_' + str(largo) + '.node', 'w')
    f.write("{} 2 0 0\n".format(largo))
    for i in range(0, len(vertices)):
        #print(i +1, vertices[i][0], vertices[i][1])
        f.write('{0} {1} {2}\n'.format(i+1, vertices[i][0],vertices[i][1]))
    f.write('\n')
    f.close()

def write_ele(num_vertices, triangles):
    largo =len(triangles)
    #for i in range(0, largo):
    #    print(i +1, triangles[i][0] +1 , triangles[i][1] +1, triangles[i][2] +1)
    f = open('input/square2x2_' + str(num_vertices) + '.ele', 'w')
    f.write("{} 3 1\n".format(largo))
    for i in range(0, largo):
        f.write( '{0} {1} {2} {3} 1\n'.format(i +1, triangles[i][0] +1 , triangles[i][1] +1, triangles[i][2] +1))
    f.write('\n')
    f.close()

if __name__ == "__main__":
    #test_quater_disk()
   # test_rectangle()
   square2x2()
