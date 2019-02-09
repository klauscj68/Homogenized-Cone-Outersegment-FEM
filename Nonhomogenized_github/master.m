function master
%Run genmesh, main_id, main_cyto
% I have two flags internal in scripts.  One is in assembla to say whether
% mass and stiffness should be built from scract.  The other is here to say
% whether the mesh should be built.

flag_build_mesh = false;
if flag_build_mesh == true
    [p_3d,f_3d,pr_3d] = genmesh;
else
    load('Mesh');
end

main_id(p_3d,f_3d);
%pause
   
main_cyto(p_3d,f_3d,pr_3d);

end

