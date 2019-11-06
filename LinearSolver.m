function x = LinearSolver( A, b )
% (Weronika Nguyen) This function performs solves a system of linear equations using LU factorization with pivoting
% and outputs x, the solution to the equation
%   This function takes a square matrix A and a vector b and outputs x such
%   that Ax=b.

    [r,c] = size(A); % r is the number of rows, c is number of columns
    % check if the matrix is square, if not then return
    if r~=c
        disp('Matrix is not square');
        return;
    end

    U=A;
    L=eye(r);
    p=(1:r)'; % this is a permutation vector

    for ii=1:c-1 %loop through cols
        % X stores the biggest value in the column, Y stores the row number
        % containing that value
            [X,Y]=max(abs(U(ii:c,ii)));
            Y=Y+ii-1; % get the right row index

            if U(Y,ii)~=0 % make sure our pivot is not equal to 0
                if Y ~= ii
                    % swap current row with row with max in U
                    curr_row = U(ii,:);
                    max_row = U(Y,:);
                    U(ii,:)=max_row;
                    U(Y,:)=curr_row;
                    
                    %swap the permutation vector
                    p([ii,Y])=p([Y ii]);
                   
                    %swap for L matrix
                    curr_l=L(ii,:);
                    max_l=L(Y,:);
                    L(ii,:)=max_l;
                    L(Y,:)=curr_l;
                   
                end
                % do the rest of LU
                for jj=(ii+1):r
                    E=eye(r);
                    E(jj,ii)=-U(jj,ii)/U(ii,ii);
                    U=E*U;
                    L=E*L;
                end
            end
    end
    % Notice we have L_inv * U = A. However L_inv isn't lower triangular.
    % We multipy both sides to get P*(L_inv)=P*A.
    % L=P*L_inv turns into a lower triangular matrix such that LU=PA.
    L = inv(L);
    L=L(p,:); % equivalent to multiplying L by P on the left side
    
    % this is a forward substitution function
    function c = forwardSub( L, b )
        [m,n] = size(L);
        c=b;
        c(1)=c(1)/L(1,1);
        for i=2:m
            for j=1:i-1
                c(i)=c(i)-L(i,j)*c(j);
            end
            c(i)=c(i)/L(i,i);
        end
    end
    
    % this is a back substitution function
    function x = backSub( U, c )
        [m,n] = size(U);
        x=c;
        x(m)=x(m)/U(m,m);
        i=m-1;
        while i>0
            j=n;
            while j>i
                x(i)=x(i)-U(i,j)*x(j);
                j=j-1;
            end
            x(i)=x(i)/U(i,i);
            i=i-1;
        end
    end
    
    c=forwardSub(L,b(p)); % compute forward substitution to calculate c
    x=backSub(U,c); % compute backward substitution and calculate x
end