function q = euler2quat(~, phi, theta, psi)
% ZYX 'frame' convention -> compact [w x y z], matching your config
q = compact(quaternion([psi theta phi], "euler", "zyx", "frame"));
end