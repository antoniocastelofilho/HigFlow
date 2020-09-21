function [NB value] = ReadDomainWriteBoundary(DIM, ND, filename, fileboundary)
    for i = 1:ND
        filenameaux = sprintf('%s-%d.amr',filename,i-1);
        file        = fopen(filenameaux,'r');
        for j = 1:DIM
           Dlow(j,i)   = fscanf(file, '%f',1);
           Dright(j,i) = fscanf(file, '%f',1);
        end
        fscanf(file, '%d',1);
        for j = 1:DIM
           DRef(j,i)   = fscanf(file, '%f',1);
        end
        fclose(file);
    end
    NB    = 0;
    eps   = 1.0e-10;
    [value Adj] = Test_Domain(DIM,ND,Dlow,Dright,eps);
    if value == 0
        display('Dominio valido');
    else
        if value == -1
            display('Dominio Invalido: blocos nao disjuntos');
            return
        else
            if value == 1
                display('Dominio Invalido: blocos disjuntos nao adjacentes');
                return
            else
                display('Dominio Invalido: blocos com adjacencia singular');
                return
            end

        end
    end
    [NB DlowBoundary DrightBoundary DRefBoundary] = Get_Boundary(ND,DIM,Dlow,Dright,DRef,Adj);
    for i = 1:NB
        filenameaux = sprintf('%s-%d.amr',fileboundary,i-1);
        file        = fopen(filenameaux,'w');
        for j = 1:DIM
            fprintf(file,'%15.10f %15.10f ',DlowBoundary(j,i),DrightBoundary(j,i));
        end
        fprintf(file,'\n');
        fprintf(file,'1\n');
        for j = 1:DIM
            fprintf(file,'%15.10f ',DRefBoundary(j,i));
        end
        fprintf(file,'1\n');
        for j = 1:DIM
            fprintf(file,'1 ');
        end
        for j = 1:DIM
            if DRefBoundary(j,i) == 0
                fprintf(file,'1 ');
            else
                fprintf(file,'%d ',round((DrightBoundary(j,i)-DlowBoundary(j,i))/DRefBoundary(j,i)));
            end
        end
        fprintf(file,'\n');
        fclose(file);
    end
    return
end


function [NB DlowBoundary DrightBoundary DRefBoundary] = Get_Boundary(n,DIM,Dlow,Dright,DRef,Adj)
    NB             = 0;
    DlowBoundary   = [];
    DrightBoundary = [];
    FaceLow        = 0;
    FaceRight      = 0;
    for j = 1:DIM
        for i = 1:n
            [NAdj FaceAdj] = Get_Face_Adj(n,i,-j,Adj);
            if NAdj == 0
                NB                   = NB+1;
                DlowBoundary(j,NB)   = Dright(j,i);
                DrightBoundary(j,NB) = Dright(j,i);
                DRefBoundary(j,NB)   = 0;
                for k = 1:j-1
                    DlowBoundary(k,NB)   = Dlow(k,i);
                    DrightBoundary(k,NB) = Dright(k,i);
                    DRefBoundary(k,NB)   = DRef(k,i);
                end
                for k = j+1:DIM
                    DlowBoundary(k,NB)   = Dlow(k,i);
                    DrightBoundary(k,NB) = Dright(k,i);
                    DRefBoundary(k,NB)   = DRef(k,i);
                end
            else
                [NGrid GridLow GridRight] = Get_Grid(n,DIM,i,-j,NAdj,FaceAdj,Dlow,Dright);
                for s = 1:NGrid
                    NB                   = NB+1;
                    DlowBoundary(j,NB)   = Dright(j,i);
                    DrightBoundary(j,NB) = Dright(j,i);
                    DRefBoundary(j,NB)   = 0;
                    for k = 1:j-1
                        DlowBoundary(k,NB)   = GridLow(k,s);
                        DrightBoundary(k,NB) = GridRight(k,s);
                        DRefBoundary(k,NB)   = DRef(k,s);
                    end
                    for k = j+1:DIM
                        DlowBoundary(k,NB)   = GridLow(k,s);
                        DrightBoundary(k,NB) = GridRight(k,s);
                        DRefBoundary(k,NB)   = DRef(k,s);
                    end
                end
            end
            [NAdj FaceAdj] = Get_Face_Adj(n,i,j,Adj);
            if NAdj == 0
                NB                   = NB+1;
                DlowBoundary(j,NB)   = Dlow(j,i);
                DrightBoundary(j,NB) = Dlow(j,i);
                DRefBoundary(j,NB)   = 0;
                for k = 1:j-1
                    DlowBoundary(k,NB)   = Dlow(k,i);
                    DrightBoundary(k,NB) = Dright(k,i);
                    DRefBoundary(k,NB)   = DRef(k,i);
                end
                for k = j+1:DIM
                    DlowBoundary(k,NB)   = Dlow(k,i);
                    DrightBoundary(k,NB) = Dright(k,i);
                    DRefBoundary(k,NB)   = DRef(k,i);
                end
            else
                [NGrid GridLow GridRight] = Get_Grid(n,DIM,i,j,NAdj,FaceAdj,Dlow,Dright);
                for s = 1:NGrid
                    NB                   = NB+1;
                    DlowBoundary(j,NB)   = Dlow(j,i);
                    DrightBoundary(j,NB) = Dlow(j,i);
                    DRefBoundary(j,NB)   = 0;
                    for k = 1:j-1
                        DlowBoundary(k,NB)   = GridLow(k,s);
                        DrightBoundary(k,NB) = GridRight(k,s);
                        DRefBoundary(k,NB)   = DRef(k,s);
                    end
                    for k = j+1:DIM
                        DlowBoundary(k,NB)   = GridLow(k,s);
                        DrightBoundary(k,NB) = GridRight(k,s);
                        DRefBoundary(k,NB)   = DRef(k,s);
                    end
                end
            end
        end
    end
    return
end


function [NAdj FaceAdj] = Get_Face_Adj (n, i, j, Adj)
    NAdj    = 0;
    FaceAdj = [];
    for k = 1:i-1
        if Adj(i,k) == j
            NAdj = NAdj + 1;
            FaceAdj = [FaceAdj k];
        end
    end
    for k = i+1:n
        if Adj(i,k) == -j
            NAdj = NAdj + 1;
            FaceAdj = [FaceAdj k];
        end
    end
    return
end

function [N GridLow GridRight] = Get_Grid(n,DIM,i,j,NAdj,FaceAdj,Dlow,Dright)
    for k = 1:DIM
        NGrid(k)     = 0;
    end
    for k = 1:abs(j)-1
        NGrid(k)          = NGrid(k) + 1;
        Grid{k}(NGrid(k)) = Dlow(k,i);
        NGrid(k)          = NGrid(k) + 1;
        Grid{k}(NGrid(k)) = Dright(k,i);
    end
    for k = abs(j)+1:DIM
        NGrid(k)          = NGrid(k) + 1;
        Grid{k}(NGrid(k)) = Dlow(k,i);
        NGrid(k)          = NGrid(k) + 1;
        Grid{k}(NGrid(k)) = Dright(k,i);
    end
    for s = 1:NAdj
        for k = 1:abs(j)-1
            if (Dlow(k,FaceAdj(s)) > Dlow(k,i)) && (Dlow(k,FaceAdj(s)) < Dright(k,i))
                NGrid(k)          = NGrid(k) + 1;
                Grid{k}(NGrid(k)) = Dlow(k,FaceAdj(s));
            end
            if (Dright(k,FaceAdj(s)) > Dlow(k,i)) && (Dright(k,FaceAdj(s)) < Dright(k,i))
                NGrid(k)          = NGrid(k) + 1;
                Grid{k}(NGrid(k)) = Dright(k,FaceAdj(s));
            end
        end
        for k = abs(j)+1:DIM
            if (Dlow(k,FaceAdj(s)) > Dlow(k,i)) && (Dlow(k,FaceAdj(s)) < Dright(k,i))
                NGrid(k)          = NGrid(k) + 1;
                Grid{k}(NGrid(k)) = Dlow(k,FaceAdj(s));
            end
            if (Dright(k,FaceAdj(s)) > Dlow(k,i)) && (Dright(k,FaceAdj(s)) < Dright(k,i))
                NGrid(k)          = NGrid(k) + 1;
                Grid{k}(NGrid(k)) = Dright(k,FaceAdj(s));
            end
        end
    end
    if j > 0
        NGrid(abs(j))               = 1;
        Grid{abs(j)}(NGrid(abs(j))) = Dlow(abs(j),i);
    else
        NGrid(abs(j))               = 1;
        Grid{abs(j)}(NGrid(abs(j))) = Dright(abs(j),i);
    end

    for k = 1:abs(j)-1
        [Grid{k}] = Reorder(NGrid(k),Grid{k});
    end
    for k = abs(j)+1:DIM
        [Grid{k}] = Reorder(NGrid(k),Grid{k});
    end
    for k = 1:abs(j)-1
        NGrid(k)  = NGrid(k)-1;
    end
    for k = abs(j)+1:DIM
        NGrid(k)  = NGrid(k)-1;
    end

    N       = 0;
    [Base]  = InitGrid(DIM,NGrid);
    for g = 0:Base(DIM+1)-1
       [FaceGrid] = GridCoords(DIM,g,Base);
       N          = N + 1;
       for k = 1:abs(j)-1
           GridLow(k,N)   = Grid{k}(FaceGrid(k)+1);
           GridRight(k,N) = Grid{k}(FaceGrid(k)+2);
       end
       for k = abs(j)+1:DIM
           GridLow(k,N)   = Grid{k}(FaceGrid(k)+1);
           GridRight(k,N) = Grid{k}(FaceGrid(k)+2);
       end
       GridLow(abs(j),N)   = Grid{abs(j)}(1);
       GridRight(abs(j),N) = Grid{abs(j)}(1);
       value               = 0;
       for s = 1:NAdj
           [in] = InFace(DIM,abs(j),GridLow(:,N),GridRight(:,N),Dlow(:,FaceAdj(s)),Dright(:,FaceAdj(s)));
           if in == 1
               value = 1;
           end
       end
       if value == 1
           N = N - 1;
       end
    end
    return
end


function [in] = InFace(DIM, j, XLow, XRight, YLow, YRight)
    in = 1;
    for k = 1:j-1
        if XLow(k) < YLow(k)
            in = 0;
            return
        end
        if XRight(k) > YRight(k)
            in = 0;
            return
        end
    end
    for k = j+1:DIM
        if XLow(k) < YLow(k)
            in = 0;
            return
        end
        if XRight(k) > YRight(k)
            in = 0;
            return
        end
    end
    return
end


function [Base] = InitGrid(n, Division)
   Base(1) = 1;
   for i = 2:n+1
      Base(i) = Base(i-1)*Division(i-1);
   end
   return
end


function [Grid] = GridCoords(n,i,Base)
   copy = i;
   for j = n:-1:2
      aux     = mod(copy,Base(j));
      Grid(j) = (copy-aux)/Base(j);
      copy    = aux;
   end
   Grid(1) = copy;
   return
end

function [x] = Reorder (n, x)
    if n > 1
        for i = 1:n-1
            for j = i+1:n
                if x(j) < x(i)
                    y    = x(j);
                    x(j) = x(i);
                    x(i) = y;
                end
            end
        end
    end
    return
end


function [value adj] = Test_Domain (DIM, n, Dlow, Dright, eps)
   for k = 1:n
      adj(k,k) = 0;
      for i = k+1:n
        for j = 1:DIM
            [lr(j) rl(j)] = Classify(Dlow(j,k),Dright(j,k),Dlow(j,i),Dright(j,i),eps);
        end
        [value] = Disjoint(DIM,lr,rl);
        if value == 1
            adj(k,i) = 0;
            adj(i,k) = 0;
        else
            [value] = Adjacent(DIM,lr,rl);
            if value ~= 0
                if value ~= DIM+1
                    adj(k,i) = value;
                    adj(i,k) = value;
                    value    = 1;
                else
                    adj(k,i) = DIM+1;
                    adj(i,k) = DIM+1;
                    value    = 2;
                end
            else
                adj(k,i) = -1;
                adj(i,k) = -1;
                value    = -1;
                return
            end
        end
      end
   end
   if value ~= 2
       [value] = One_CC(n,adj);
   end
   return
end


function [value] = One_CC (n, adj)
    comp    = zeros(n,1);
    comp(1) = 1;
    for k = 1:n-1
        for i = k+1:n
            if (comp(k) == 1) && (comp(i) == 0)
                if adj(k,i) ~= 0
                    comp(i) = 1;
                end
            end
        end
    end
    value = 0;
    for i = 1:n
        if comp(i) == 0
            value = 1;
        end
    end

    return
end


function [value] = Disjoint (DIM, x, y)
    value = 1;
    for j = 1:DIM
        if x(j) ==  1
            return
        end
        if y(j) == -1
            return
        end
    end
    value = 0;
    return
end


function [value] = Adjacent (DIM, x, y)
    for j = 1:DIM
        if x(j)*y(j) == 0
            if DIM > 1
                XX = [x(1:j-1) x(j+1:DIM)];
                YY = [y(1:j-1) y(j+1:DIM)];
                [value] = Disjoint(DIM-1,XX,YY);
                aux = 1;
                for i = 1:DIM-1
                    aux = XX(i)*YY(i);
                end
                if aux ~= 0
                    if value == 0
                        if x(j) == 0
                            value = -j;
                            return
                        else
                            value = j;
                            return
                        end
                    end
                else
                    value = DIM+1;
                    return
                end
            else
                if x(j) == 0
                    value = -j;
                    return
                else
                    value = j;
                    return
                end
            end
        end
    end
    value = 0;
    return
end


function [lr rl] = Classify (xl, xr, yl, yr, eps)
    [lr] = Test_Point(xl,yr,eps);
    [rl] = Test_Point(xr,yl,eps);
    return
end


function [value] = Test_Point (x, y, eps)
    if abs(x-y) < eps
        value = 0;
    else
        if x < y
            value = -1;
        else
            value = 1;
        end
    end
    return
end