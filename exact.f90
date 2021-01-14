subroutine exact(Nx, Ny, rho_exact, Gauss, x, y, time)
    use global
    implicit none
    integer :: Nx, Ny
    real :: local_x(k), x(Nx+1), local_y(k), y(Ny+1), dx, dy
    real :: Gauss(k), rho_exact(dim*k*k,Nx*Ny), time
    integer :: i, j, ii, jj
    real, external :: fun, dfun
    real :: u, f, df

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
                ! 线性对流方程精确解
                rho_exact(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 100.0D0 + SIN(local_x(ii)-pi*time) + SIN(local_y(jj)-pi*time)
                rho_exact(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 100.0D0 + pi*(SIN(local_x(ii)-pi*time) + SIN(local_y(jj)-pi*time))
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
                ! 线性对流方程精确解
                rho_exact(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 10.0D0 + exp(-time)*(SIN(local_x(ii)) + SIN(local_y(jj)))
                rho_exact(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 10.0D0 - exp(-time)*(COS(local_x(ii)) + COS(local_y(jj)))
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
                !rho_exact(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = -2.0D0 - SIN(local_x(ii)+local_y(jj)-2.0*time)
                !rho_exact(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = -2.0D0 - SIN(local_x(ii)+local_y(jj)-2.0*time)
                
                rho_exact(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 10.0D0 + SIN(local_x(ii)-pi*time) + SIN(local_y(jj)-pi*time)
                rho_exact(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 10.0D0 + pi*(SIN(local_x(ii)-pi*time) + SIN(local_y(jj)-pi*time))
            end do
            end do
        end do
        end do      
        
    end select

    return
end subroutine exact
    
function fun(u,x,y,t)
    !Burgers 方程精确解Newton迭代 非线性方程
    implicit none
    real :: u, x, y, t, pi = DACOS(-1.0D0), fun
    fun = u - sin(2.0*pi*(x+y-2.0*u*t))/2.0D0
    return
end function
    
function dfun(u,x,y,t)
    !Burgers 方程精确解Newton迭代 非线性方程偏导
    implicit none
    real :: u, x, y, t, pi = DACOS(-1.0D0), dfun
    dfun = 1.0 + 2.0D0*pi*t*cos(2.0*pi*(x+y-2.0*u*t))
    return
end function