

######
## Element class tags

fourNodeEleTags = [31, 32, 52, 53, 55, 156, 157, 173, 174, 203]
triNodeEleTags = [33, 167, 168, 204]
eightNodeBrickEleTags = [36, 38, 43, 44, 56]
MVLEMEleTags = [162, 164, 212, 213]
tetEleTags = [179]
twoNodeLink = [86]


# Add beam column elements
#
twoNodeBeamColumn = [[i for i in range(1,31)]+[i for i in range(62,80)]+[34,35,621,63,64,640,641,642,731,751,30766,30765]][0]


####
### Element tag groupf for stress and strain recorder
####
shell_4N4GP_EleTags = [31, 32, 52, 53, 55, 156, 157, 173, 174, 203, 212, 213]  # Quad only forks for 2d/ fix assignment?
shell_3N3GP_EleTags = [33, 167, 168, 204]
tetra_4N4GP_EleTags = [179]
brick_8N8GP_EleTags = [36, 38, 56]  # Bricks only had 3 dofs; do not mix with others










