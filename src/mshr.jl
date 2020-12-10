############################################
#
#    mshr
#
# https://bitbucket.org/dolfin-project/mshr/wiki/Home
#
############################################

# no need to export
# not exported in `mshr`
@pyclass mshr CSGGeometry

# 2D primitives
@pyclass mshr Rectangle CSGGeometry
@pyclass mshr Circle CSGGeometry
@pyclass mshr Polygon CSGGeometry
@pyclass mshr Ellipse CSGGeometry

# 3D primitives
@pyclass mshr Cylinder CSGGeometry
@pyclass mshr Box CSGGeometry
@pyclass mshr Surface3D CSGGeometry
@pyclass mshr Cone CSGGeometry
@pyclass mshr Ellipsoid CSGGeometry
@pyclass mshr Sphere CSGGeometry
@pyclass mshr Tetrahedron CSGGeometry

@pyfunc mshr generate_mesh

+(geo1::CSGGeometry, geo2::CSGGeometry) = CSGGeometry(geo1.pyobject + geo2.pyobject)
-(geo1::CSGGeometry, geo2::CSGGeometry) = CSGGeometry(geo1.pyobject - geo2.pyobject) 

export Rectangle, Circle, Polygon, Ellipse,
       Cylinder, Box, Surface3D, Cone, Ellipsoid, Sphere, Tetrahedron, 
       generate_mesh



