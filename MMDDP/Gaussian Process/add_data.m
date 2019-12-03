xo1=[3.14,3.15];
xo2=[3.08,2.98];
xo3=[0.4,0.1,-0.2];
xo4=[-0.8,-2, -4];
uo1=[-26.5,-26,-25];
uo2=[-33,-28];
for a=1:2
    for b=1:2
        for c=1:3
            for d=1:3
                for e=1:3
                    for f=1:2
                        x_norm = [x_norm;[xo1(a),xo2(b),xo3(c),xo4(d)]];
                        u_norm = [u_norm;[uo1(e),uo2(f)]];
                    end
                end
            end
        end
    end
end

xo1=[0.15,-0.2];
xo2=[0.15,-0.2];
xo3=[0.15,-0.2];
xo4=[0.15,-0.2];
uo1=[2.7,3,3.4];
uo2=[2.7,3,3.4];
for a=1:2
    for b=1:2
        for c=1:2
            for d=1:2
                for e=1:3
                    for f=1:3
                        x_norm = [x_norm;[xo1(a),xo2(b),xo3(c),xo4(d)]];
                        u_norm = [u_norm;[uo1(e),uo2(f)]];
                    end
                end
            end
        end
    end
end
