%function compute_spectrum()
alpha=[-1,-0.5,0,0.5];
geo="emboite";
ni = 3;
type = 1;

matrix_P="P_"+geo+"_"+ni+"_"+type;
matrix_MTF=strings(1,size(alpha,2)-1);
load("../output/matrices/P_"+geo+"_"+ni+"_"+type+".mat");
nrows = size(eval(matrix_P),1);
matrix = zeros(nrows,size(alpha,2)*2);
headers = [];
    
    for i=0:size(alpha,2)-1
        disp(ni+" "+geo+" "+type+" "+alpha(i+1))
        load("../output/matrices/MTF_"+geo+"_"+ni+"_"+type+"_"+i+".mat");
        matrix_MTF(i+1)="MTF_"+geo+"_"+ni+"_"+type+"_"+i;
        [V,D]=eig(eval(matrix_MTF(i+1)),eval(matrix_P));
                
        z_real = real(diag(D));
        z_imag = imag(diag(D));
        matrix(:,i*2+1)=z_real(:);
        matrix(:,(i+1)*2)=z_imag(:);
        headers = [headers [alpha(i+1)+"_Real",alpha(i+1)+"_Imag"]];
        
    end
    
headers = strjoin(headers, ',');

%write header to file
fid = fopen("../output/csv/spectrum_"+geo+"_"+ni+"_"+type+".csv",'w'); 
fprintf(fid,'%s\n',headers);
fclose(fid);

%write data to end of file
dlmwrite("../output/csv/spectrum_"+geo+"_"+ni+"_"+type+".csv",matrix,'-append');

