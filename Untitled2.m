shap5File = 'C:\Users\pete4\Downloads\CG_SPC_SHAP5_012K_CART.WRL';

shap5Verts = dlmread(shap5File, '\t', [6,0,20,0]);
%%
fID = fopen(shap5File);
startRow = 7;
endRow = 5867;
nRows = endRow - startRow;
Verts = textscan(fID, '%f', 'HeaderLines',6);
Verts = reshape(Verts{1}, 3, []);
Faces = textscan(fID, '%d', 'HeaderLines',3);
Faces = reshape(Faces{1}, 4, []);
Faces = double(Faces(1:3,:)) +1 ;
fclose(fID);

%%
figure;
patch('Faces', Faces', 'Vertices', Verts', ...
    'EdgeColor', 'none','FaceColor',[0.5,0.5,0.5]);
%%
rVert = sqrt(sum(Verts.^2, 1));
std(rVert)
std(rVert)/mean(rVert)