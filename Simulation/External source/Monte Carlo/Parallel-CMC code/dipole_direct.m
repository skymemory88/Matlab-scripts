function D=dipole_direct(rij,L,N)
% D=dipole_direct(rij,L,N)
% Calculates the dipole sum in a brute force way over N(i) replicas of the
% box specified by L in the i direction. This include the Lorentz factor
% thus the resulting matrix is an approximation of the infinite sum. rij
% must be a 3x1 vector, L a 3x1 vector and N a 3x1 vector of integers.
% rij and L are of units of length. The dipole matrix D is in units of
% inverse volume.