# sixDoFRigidBodyMotion_separatedLateralMotion

Repository to manage new feature in the OpenFOAM/v2406 sixDoFRigidBodyMotion
library enabling different inner and outer distance for cell stretching and
compression in lateral direction, in vertical direction and for rotations.

Working title of modified library: sixDoFRigidBodyModifiedMeshMotion

For parameters required to employ modification, see the test case
test/moveDynamicMesh/modified/constant/dynamicMeshDict.

Verified with test/floatingObject/{original,modified} that usage of modified
library does not alter cases where new parameters are not included in
dynamicMeshDict.

Developed by Johannes Palm and Claes Eskilsson.

Contributors: Johan Roenby, Sithik Aliyar, Sarath Karumathil, Amirhossein Taran