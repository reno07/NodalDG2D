subroutine initial(Nx, Ny, rho0, Gauss, x, y)
    use global
    implicit none
    integer :: Nx, Ny
    real :: Gauss(k), local_x(k), x(Nx+1), local_y(k), y(Ny+1)
    real :: rho0(dim*k*k,Nx*Ny)
    integer :: i, j, ii, jj
    real :: dx, dy

    select case (problem)
    case (1)
        do j = 1,Ny
        do i = 1,Nx
        dx = x(i+1) - x(i)
        dy = y(j+1) - y(j)
        local_x = 0.5*(x(i) + x(i+1)) + 0.5*dx*Gauss
        local_y = 0.5*(y(j) + y(j+1)) + 0.5*dy*Gauss
            do jj = 1,k
            do ii = 1,k
                ! 线性对流方程初值
                rho0(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 100.0D0 + SIN(local_x(ii)) + SIN(local_y(jj))     
                rho0(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 100.0D0 + pi*(SIN(local_x(ii)) + SIN(local_y(jj)))      
            end do
            end do
        end do
        end do
        
    case (2)
        do j = 1,Ny
        do i = 1,Nx
        dx = x(i+1) - x(i)
        dy = y(j+1) - y(j)
        local_x = 0.5*(x(i) + x(i+1)) + 0.5*dx*Gauss
        local_y = 0.5*(y(j) + y(j+1)) + 0.5*dy*Gauss
            do jj = 1,k
            do ii = 1,k
                ! 线性对流方程初值
                rho0(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 10.0D0 + SIN(local_x(ii)) + SIN(local_y(jj))     
                rho0(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 10.0D0 - COS(local_x(ii)) - COS(local_y(jj))     
            end do
            end do
        end do
        end do
        
    case (3)
        do j = 1,Ny
        do i = 1,Nx
            dx = x(i+1) - x(i)
            dy = y(j+1) - y(j)
            local_x = 0.5*(x(i) + x(i+1)) + 0.5*dx*Gauss
            local_y = 0.5*(y(j) + y(j+1)) + 0.5*dy*Gauss
            do jj = 1,k
            do ii = 1,k
                ! 带源项的非线性方程精确解
                !rho0(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = -2.0D0 - SIN(local_x(ii)+local_y(jj))
                !rho0(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = -2.0D0 - SIN(local_x(ii)+local_y(jj))
                
                rho0(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 10.0D0 + SIN(local_x(ii)) + SIN(local_y(jj))     
                rho0(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 10.0D0 + pi*(SIN(local_x(ii)) + SIN(local_y(jj)))  
            end do
            end do
        end do
        end do 
        
    case (4)
        do j = 1,Ny
        do i = 1,Nx
            dx = x(i+1) - x(i)
            dy = y(j+1) - y(j)
            local_x = 0.5*(x(i) + x(i+1)) + 0.5*dx*Gauss
            local_y = 0.5*(y(j) + y(j+1)) + 0.5*dy*Gauss
            do jj = 1,k
            do ii = 1,k
                ! 线性对流方程精确解
                rho0(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 0.5D0
                rho0(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 0.5D0*(1.0-tanh((sqrt(local_x(ii)**2+local_y(jj)**2)/2.0D0-0.5D0)/0.2D0))
            end do
            end do
        end do
        end do 
        
    end select


    return
end subroutine initial