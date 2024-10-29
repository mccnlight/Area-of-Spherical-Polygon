import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

def spherical_to_cartesian(r, theta, phi):
    """
    Convert spherical coordinates (r, theta, phi) to Cartesian coordinates (x, y, z).
    """
    x = r * math.sin(phi) * math.cos(theta)  # Calculate x coordinate
    y = r * math.sin(phi) * math.sin(theta)  # Calculate y coordinate
    z = r * math.cos(phi)                     # Calculate z coordinate
    coordinates = [x, y, z]

    # Adjust very small values close to zero to exactly zero
    for i in range(3):
        if coordinates[i] < 0.00000001 and coordinates[i] > 0:
            coordinates[i] = 0
    return coordinates

def dot_product(v1, v2):
    """
    Calculate the dot product of two vectors v1 and v2.
    """
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

def magnitude(v):
    """
    Calculate the magnitude (length) of a vector v.
    """
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def angle_between(v1, v2):
    """
    Calculate the angle in radians between two vectors v1 and v2 using the dot product.
    """
    dot_prod = dot_product(v1, v2)             # Calculate dot product
    mag_v1 = magnitude(v1)                      # Magnitude of vector v1
    mag_v2 = magnitude(v2)                      # Magnitude of vector v2
    cos_theta = dot_prod / (mag_v1 * mag_v2)  # Calculate cosine of the angle
    angle = math.acos(cos_theta)               # Calculate angle in radians
    if angle != 0:
        return angle
    else:
        return math.radians(180)                # Return 180 degrees if the angle is zero

def points_to_normal(p1, p2, p3):
    """
    Compute the normal vector to the plane formed by three points p1, p2, p3 in 3D space.
    """
    x1, y1, z1, x2, y2, z2, x3, y3, z3 = *p1, *p2, *p3  # Unpack points
    a = (x2 - x1, y2 - y1, z2 - z1)                    # Vector from p1 to p2
    b = (x3 - x1, y3 - y1, z3 - z1)                    # Vector from p1 to p3
    # Calculate normal vector using cross product
    normal = (
        (a[1] * b[2]) - (b[1] * a[2]),
        -1 * ((a[0] * b[2]) - (b[0] * a[2])),
        (a[0] * b[1]) - (b[0] * a[1])
    )
    return normal

def spherical_polygon_angles(vertices):
    """
    Calculate the interior angles of a polygon defined by vertices in spherical coordinates.
    """
    vertices = list(set(vertices))  # Remove duplicate vertices
    vertices.sort(key=lambda v: (v[1], v[2]))  # Sort vertices by theta and phi
    n = len(vertices)
    angles = []

    for i in range(n):
        origin = (0, 0, 0)  # Define the origin point
        # Convert spherical coordinates to Cartesian for each vertex
        p1 = spherical_to_cartesian(*vertices[i])
        p2 = spherical_to_cartesian(*vertices[(i + 1) % n])
        p3 = spherical_to_cartesian(*vertices[(i + 2) % n])

        vector1 = points_to_normal(origin, p1, p2)  # Get normal vector for triangle formed by p1 and p2
        vector2 = points_to_normal(origin, p2, p3)  # Get normal vector for triangle formed by p2 and p3

        angle = angle_between(vector1, vector2)  # Calculate angle between two normal vectors
        angles.append(angle)  # Store the angle
    return angles

def area_of_the_polygon(interior_angles, radius):
    """
    Calculate the area of a polygon on a sphere given its interior angles and the radius of the sphere.
    """
    sum_of_the_angles = sum(interior_angles)  # Sum of the interior angles
    # Calculate the excess of the angles over (n-2)*π
    excess = sum_of_the_angles - (len(interior_angles) - 2) * math.pi
    area = (radius**2) * excess  # Area of the polygon
    return area

def visualize_sphere_and_polygon(vertices, radius, elev=30, azim=30):
    """
    Visualize the sphere and the polygon defined by its vertices in 3D space.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Create data for a sphere
    u = np.linspace(0, 2 * np.pi, 100)  # Azimuthal angle
    v = np.linspace(0, np.pi, 100)      # Polar angle
    x = radius * np.outer(np.cos(u), np.sin(v))  # X coordinates
    y = radius * np.outer(np.sin(u), np.sin(v))  # Y coordinates
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))  # Z coordinates

    # Plot the sphere
    ax.plot_surface(x, y, z, color='c', alpha=0.3, rstride=4, cstride=4)

    # Plot the polygon vertices
    polygon_points = []
    for vertex in vertices:
        cartesian = spherical_to_cartesian(*vertex)  # Convert each vertex to Cartesian
        polygon_points.append(cartesian)
        ax.scatter(cartesian[0], cartesian[1], cartesian[2], color='r', s=50)  # Plot vertex points

    # Close the polygon by connecting the last point to the first
    polygon_points.append(polygon_points[0])

    # Convert list to numpy array for easier plotting
    polygon_points = np.array(polygon_points)

    # Plot the polygon surface
    poly3d = [[polygon_points[i] for i in range(len(polygon_points))]]
    ax.add_collection3d(Poly3DCollection(poly3d, facecolors='blue', linewidths=1, edgecolors='r', alpha=.25))

    # Plot the polygon edges
    ax.plot(polygon_points[:, 0], polygon_points[:, 1], polygon_points[:, 2], color='b', lw=2)

    # Plot axes
    axes_length = radius * 1.5
    ax.plot([0, axes_length], [0, 0], [0, 0], color='k', linewidth=2)  # X-axis
    ax.plot([0, 0], [0, axes_length], [0, 0], color='k', linewidth=2)  # Y-axis
    ax.plot([0, 0], [0, 0], [0, axes_length], color='k', linewidth=2)  # Z-axis

    # Set the camera view
    ax.view_init(elev=elev, azim=azim)

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Polygon on a Sphere')
    plt.show()  # Display the plot

# Define the radius of the sphere
radius = 10

# Define the polygon vertices in spherical coordinates (r, theta, phi)
polygon_vertices = [
    (radius, math.radians(0), math.radians(0)),
    (radius, math.radians(0), math.radians(90)),
    (radius, math.radians(90), math.radians(90))
]

# Ensure there are enough points to form a polygon
if len(polygon_vertices) < 3:
    print("At least 3 points are needed to form a polygon!")
else:
    # Calculate and print the area of the polygon
    print('Area of the Polygon: ', area_of_the_polygon(spherical_polygon_angles(polygon_vertices), radius))
    # Visualize the sphere and the polygon
    visualize_sphere_and_polygon(polygon_vertices, radius)
