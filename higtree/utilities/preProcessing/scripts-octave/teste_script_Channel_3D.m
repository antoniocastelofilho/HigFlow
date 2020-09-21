DIM = 3;

filedomain   = 'Channel_3d_d/mesh-channel3D-d';
fileboundary = 'Channel_3d_bc/mesh-channel3D-bc';
fileprobe    = 'Channel_3d_probe/mesh-channel3D-probe';
filescript   = 'Channel_3d_script/mesh-channel3D-script';
ND           = 16;
NP           = 0;
delta        = 0.01;
reynolds     = 10.0;
[NB value]   = ReadDomainWriteBoundary(DIM,ND,filedomain,fileboundary);
if value == 0
    WriteScript(DIM,ND,NB,NP,delta,reynolds,filedomain,fileboundary,fileprobe,filescript);
end