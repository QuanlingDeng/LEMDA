%% RMSE root mean square error and PCC pattern cross correlation
rmse = zeros(ndim,1); rrmse = zeros(ndim,1); pcc = zeros(ndim,1); 

nr = length(daSol(1,:));  nl = round(nr/2); %nn=nn/2;
np = length(nl:nr);

for j=0:5
    vec1 = trueData(j*Lo+1:j*Lo+L,nl:nr);
    vec2 = daSol(j*L+1:j*L+L,nl:nr); 
    
    if j<2
        vec2 = (abs(vec1 - vec2) <= pi) .* vec2 + (vec1 - vec2 > pi) .* vec2 + (vec2 - vec1 > pi) .* (vec2 - 2*pi);
        vec1 = (abs(vec1 - vec2) <= pi) .* vec1 + (vec1 - vec2 > pi) .* (vec1 - 2*pi) + (vec2 - vec1 > pi) .* vec1;
    end
    rmse(j*L+1:j*L+L) = sum( (vec1 - vec2).^2, 2)/np;
    rrmse(j*L+1:j*L+L) = sum( (vec1).^2, 2)/np;
    pcc(j*L+1:j*L+L) = dot(vec1,vec2,2) ./ sqrt(dot(vec1,vec1,2) .* dot(vec2,vec2,2) );
end

vec1 = trueData(6*Lo+1:6*Lo+Dim_UB,nl:nr);
vec2 = daSol(6*L+1:6*L+Dim_UB,nl:nr);
rmse(6*L+1:end) = sum( (vec1 - vec2).*conj(vec1 - vec2), 2)/np;
rrmse(6*L+1:end) = sum( (vec1).*conj(vec1), 2)/np;
pcc(6*L+1:end) = dot(vec1,vec2,2) ./ sqrt(dot(vec1,vec1,2) .* dot(vec2,vec2,2) );
pcc = real(pcc); rrmse = real(rrmse);

rrmse = sqrt(rmse./rrmse); rmse = sqrt(rmse); 

err = [mean(rrmse(1:L)) mean(rrmse(L+1:2*L)) mean(rrmse(2*L+1:3*L)) mean(rrmse(3*L+1:4*L)) mean(rrmse(4*L+1:5*L)) mean(rrmse(5*L+1:6*L)) mean(rrmse(6*L+1:end))];
cc = [mean(pcc(1:L)) mean(pcc(L+1:2*L)) mean(pcc(2*L+1:3*L)) mean(pcc(3*L+1:4*L)) mean(pcc(4*L+1:5*L)) mean(pcc(5*L+1:6*L)) mean(pcc(6*L+1:end))];
[err; cc]

%save(['./da/SupF', num2str(Lo), 'TO', num2str(Ll),'S',num2str(Ls),'NS',num2str(round(10*param.sigOcn)),'T',num2str(T),'.mat'],'err','cc','obsData', 'daSol');%'trueT', 'daSV',

