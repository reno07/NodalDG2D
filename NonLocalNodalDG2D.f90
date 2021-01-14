!**************************************************************************
! 2D + nodal LDG + SSP-RK
! 求解 非局域热传导模型
!**************************************************************************
module global
    implicit none
    real, parameter :: pi = DACOS(-1.0D0)
    integer, parameter :: MAXITER = 100000  ! 最大迭代次数
    integer, parameter :: dim = 2
    integer, parameter :: k = 4     ! GaussLobatto积分点数量 得到 (k-1) 元
    integer, parameter :: kt = 2    ! SSP Runge-Kutta 阶数 
    integer, parameter :: Loop = 3
    integer, parameter :: problem = 3 ! 1-线性对流 2-线性热传导方程 3-SKT精度测试
    real :: xl, xr, yb, yu, T, lambda, para_k
    integer, allocatable :: n1(:),n2(:),BoundaryIndex(:),XFluxIndex(:),YFluxIndex(:)
    
    contains
    subroutine CreateFluxPara
    allocate(BoundaryIndex(4*k-4), XFluxIndex(4*k-4), YFluxIndex(4*k-4), n1(4*k-4), n2(4*k-4))
        select case (k)
        case (3)
            BoundaryIndex = (/1,2,3,4,6,7,8,9/)
            XFluxIndex = (/3,5,1,6,4,9,5,7/)
            YFluxIndex = (/7,8,9,5,5,1,2,3/)
            n1 = (/-1,0,1,-1,1,-1,0,1/)  
            n2 = (/-1,-1,-1,0,0,1,1,1/)
        case (4)
            BoundaryIndex = (/1,2,3,4,5,8,9,12,13,14,15,16/)
            XFluxIndex = (/4,1,1,1,8,5,12,9,16,1,1,13/)
            YFluxIndex = (/13,14,15,16,1,1,1,1,1,2,3,4/)
            n1 = (/-1,0,0,1,-1,1,-1,1,-1,0,0,1/)  
            n2 = (/-1,-1,-1,-1,0,0,0,0,1,1,1,1/)
        case (5)
            BoundaryIndex = (/1,2,3,4,5,6,10,11,15,16,20,21,22,23,24,25/)
            XFluxIndex = (/5,1,1,1,1,10,6,15,11,20,16,25,1,1,1,21/)
            YFluxIndex = (/21,22,23,24,25,1,1,1,1,1,1,1,2,3,4,5/)
            n1 = (/-1,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,1/)
            n2 = (/-1,-1,-1,-1,-1,0,0,0,0,0,0,1,1,1,1,1/)
        end select
    end subroutine
    
    subroutine DeleteFluxPara
        deallocate(BoundaryIndex, XFluxIndex, YFluxIndex, n1, n2)
    end subroutine
    
    subroutine CreateDomainPara
        select case (problem)
        case (1)
            xl = 0.D0       ! 计算区间左端点
            xr = 2.0*pi     ! 计算区间右端点
            yb = 0.D0       ! 计算区间下端点
            yu = 2.0*pi     ! 计算区间上端点
            T = 0.1         ! 输出时间
            lambda = 0.0
            para_k = pi*pi
        case (2)
            xl = 0.D0       ! 计算区间左端点
            xr = 2.0*pi     ! 计算区间右端点
            yb = 0.D0       ! 计算区间下端点
            yu = 2.0*pi     ! 计算区间上端点
            T = 0.1         ! 输出时间
            lambda = 2.0
            para_k = 1.0
        case (3)
            xl = 0.D0       ! 计算区间左端点
            xr = 2.0*pi     ! 计算区间右端点
            yb = 0.D0       ! 计算区间下端点
            yu = 2.0*pi     ! 计算区间上端点
            T = 0.05        ! 输出时间
            lambda = 0.0
            para_k = 0.0
        end select
    end subroutine

end module global
    
subroutine PartialLagrangian(Gauss, k, D)
    implicit none
    integer :: k
    real :: Gauss(k), D(k,k)
    integer :: i,j,l,ii
    real :: denominator, numerator, sub_numerator
    do j = 1,k
        do l = 1,k
            denominator = 1.0
            numerator = 0.0
            do i = 1,k
                if (i /= l) then
                    denominator = denominator * (Gauss(l) - Gauss(i))
                    sub_numerator = 1
                    do ii = 1,k
                        if (ii/=i .and. ii/=l) then
                            sub_numerator = sub_numerator * (Gauss(j) - Gauss(ii))
                        end if
                    end do
                    numerator = numerator + sub_numerator
                end if
            end do
            D(j,l) = numerator / denominator
        end do
    end do
    return
    
end subroutine PartialLagrangian
    
subroutine Kronecker(p1,q1,A,p2,q2,B,AB)
! 输入m1xn1的矩阵A m2xn2的矩阵B 输出矩阵AB
    use global
    implicit none
    integer :: i, j, p1, q1, p2, q2
    real :: A(p1,q1), B(p2,q2), AB(p1*p2, q1*q2)
    
    do i = 1,p1
        do j = 1,q1
            AB((i-1)*p2+1:i*p2, (j-1)*q2+1:j*q2) = A(i,j)*B
        end do
    end do
    return
end subroutine Kronecker

subroutine GaussLobatto_1d(k,Gauss,wt)
    implicit none
    integer :: k
    real :: Gauss(k), wt(k)
    select case (k)
    case (3)
        Gauss = (/ -1.D0, 0.D0, 1.D0 /)
        wt = (/ 1.D0/3.D0, 4.D0/3.D0, 1.D0/3.D0 /)
    case (4)
        Gauss = (/ -1.D0, -sqrt(1.D0/5.D0), sqrt(1.D0/5.D0), 1.D0 /)
        wt = (/ 1.D0/6.D0, 5.D0/6.D0, 5.D0/6.D0, 1.D0/6.D0 /)
    case (5)
        Gauss = (/ -1.0D0,-DSQRT(3.0D0/7.0D0), 0.0D0, DSQRT(3.0D0/7.0D0), 1.0D0 /)
        wt = (/ 1.0D0/10.0D0, 49.0D0/90.0D0, 32.0D0/45.0D0, 49.0D0/90.0D0, 1.0D0/10.0D0 /)
    end select
    return
    end subroutine
    
subroutine SSP_RK(Nx, Ny, rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs, dt, time)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(k*k*dim, dim*k*k) :: BigMat1, BigMat2, BigMat3, BigMat4
    real :: rho(dim*k*k,Nx*Ny), rhs(dim*k*k,Nx*Ny), rho1(dim*k*k,Nx*Ny), rho2(dim*k*k,Nx*Ny), s(dim*k*k,Nx*Ny)
    real :: x(Nx+1), y(Ny+1), Gauss(k), wt(k), dt, time
    
    select case (kt)
    case (1)
        call SpaceDis(Nx, Ny, rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        rho = rho + dt*rhs
    case (2)
        !call SpaceDis(Nx, Ny, rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        !rho1 = rho + dt*rhs
        !call SpaceDis(Nx, Ny, rho1, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        !rho = 0.5D0*rho + 0.5D0*(rho1 + dt*rhs)
        
        call func_s(Nx, Ny, Gauss, x, y, time, rho, s)
        call SpaceDis(Nx, Ny, rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        rho1 = rho + dt*(rhs+s)
        call func_s(Nx, Ny, Gauss, x, y, time+dt, rho1, s)
        call SpaceDis(Nx, Ny, rho1, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)        
        rho = 0.5D0*rho + 0.5D0*(rho1 + dt*(rhs+s))
    case (3)
        !call SpaceDis(Nx, Ny, rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        !rho1 = rho + dt*rhs
        !call SpaceDis(Nx, Ny, rho1, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        !rho2 = 0.75*rho + 0.25*(rho1 + dt*rhs)
        !call SpaceDis(Nx, Ny, rho2, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        !rho = 1.0D0/3.0D0*rho + 2.0D0/3.0D0*(rho2 + dt*rhs)
        
        
        call SpaceDis(Nx, Ny, rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        call func_s(Nx, Ny, Gauss, x, y, time, rho, s)
        rho1 = rho + dt*(rhs+s)
        call SpaceDis(Nx, Ny, rho1, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        call func_s(Nx, Ny, Gauss, x, y, time+2.0*dt/3.0, rho1, s)
        rho2 = 0.75*rho + 0.25*(rho1 + dt*(rhs+s))
        call SpaceDis(Nx, Ny, rho2, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs)
        call func_s(Nx, Ny, Gauss, x, y, time+dt, rho2, s)
        rho = 1.0D0/3.0D0*rho + 2.0D0/3.0D0*(rho2 + dt*(rhs+s))
    end select
    return
end subroutine SSP_RK
    
subroutine func_s(Nx, Ny, Gauss, x, y, time, rho, s)
    use global
    implicit none
    integer :: Nx, Ny
    real :: local_x(k), x(Nx+1), local_y(k), y(Ny+1), dx, dy
    real :: Gauss(k), s(dim*k*k,Nx*Ny), time, A, B, rho(dim*k*k,Nx*Ny)
    integer :: i, j, ii, jj
    
    do j = 1,Ny
    do i = 1,Nx
        dx = x(i+1) - x(i)
        dy = y(j+1) - y(j)
        local_x = 0.5*(x(i) + x(i+1)) + 0.5*dx*Gauss
        local_y = 0.5*(y(j) + y(j+1)) + 0.5*dy*Gauss
        do jj = 1,k
        do ii = 1,k
            ! SKT 精度测试源项
            !A = local_x(ii)+local_y(jj)-2.0*time
            !B = local_x(ii)+local_y(jj)-2.0*time
            !s(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 0.0
            !s(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 2.0*COS(B)-2.0*lambda*SIN(B)-2.0*COS(A)*SIN(A)
            
            
            A = local_x(ii)-pi*time
            B = local_y(jj)-pi*time
            s(((jj-1)*k+ii)*dim-1, (j-1)*Nx+i) = 0.0
            !s(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = lambda*pi*(SIN(A)+SIN(B)) + (COS(A)+COS(B))*((SIN(A)+SIN(B))**1+20.0-pi**2)
            s(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = lambda*pi*(SIN(A)+SIN(B)) + (COS(A)+COS(B))*(rho(((jj-1)*k+ii)*dim-1,(j-1)*Nx+i)-pi**2)
            !s(((jj-1)*k+ii)*dim, (j-1)*Nx+i) = 0.0
        end do
        end do
        !write(*,*) "I am here!"
    end do
    end do  
    return
end subroutine func_s

subroutine func_f1(Nx, Ny, f1, rho)
    use global
    implicit none
    integer :: Nx, Ny
    real :: f1(dim*k*k,Nx*Ny), rho(dim*k*k,Nx*Ny)
    ! 线性对流方程
    f1 = rho
    !Burgers 方程
    !f1 = rho**2/2.0D0
    ! Burgers 方程 间断测试
    !f1 = rho**2
    
    return
end subroutine func_f1
    
subroutine func_f2(Nx, Ny, f2, rho)
    use global
    implicit none
    integer :: Nx, Ny
    real :: f2(dim*k*k,Nx*Ny), rho(dim*k*k,Nx*Ny)
    ! 线性对流方程
    f2 = rho
    !Burgers 方程
    !f2 = rho**2/2.0D0
    ! Burgers 方程 间断测试
    !f2 = rho**2
    
    return
end subroutine func_f2
    
subroutine func_Q(Nx, Ny, theta, rho, Q)
    use global
    implicit none
    integer :: Nx, Ny, j
    real :: Q(dim*k*k,Nx*Ny), theta(dim*k*k,Nx*Ny), rho(dim*k*k,Nx*Ny)
    ! 线性对流方程
    !do j = 1,k*k
    !    Q(dim*j-1,:) = rho(dim*j,:)
    !    Q(dim*j,:) = pi*pi*rho(dim*j-1,:)
    !enddo
    
    ! 线性对流-扩散方程
    !do j = 1,k*k
    !    Q(dim*j-1,:) = -rho(dim*j,:)
    !    Q(dim*j,:) = lambda*theta(dim*j,:) - para_k*rho(dim*j-1,:)
    !enddo
    
    ! 带源项非线性对流扩散方程
    do j = 1,k*k
        Q(dim*j-1,:) = -rho(dim*j,:)
        Q(dim*j,:) = lambda*theta(dim*j,:) - rho(dim*j-1,:)**2/2.0
    enddo

    return
end subroutine func_Q  
    
subroutine func_G(Nx, Ny, rho, G)
    use global
    implicit none
    integer :: Nx, Ny, j
    real :: G(dim*k*k,Nx*Ny), rho(dim*k*k,Nx*Ny)
    ! 线性对流方程
    do j = 1,k*k
        G(dim*j-1,:) = rho(dim*j,:)
        G(dim*j,:) = pi*pi*rho(dim*j-1,:)
    enddo

    return
end subroutine func_G   

function periodic_boundary(i,N)
    implicit none
    integer :: i, N, periodic_boundary
    periodic_boundary = i
    
    if (i < 1) then
        periodic_boundary = N
    else if (i > N) then
        periodic_boundary = 1
    else
        periodic_boundary = i
    end if
    return
end function periodic_boundary
    
function riemann_boundary(i,N)
    implicit none
    integer :: i, N, riemann_boundary
    riemann_boundary = i
    
    if (i < 1) then
        riemann_boundary = 1
    else if (i > N) then
        riemann_boundary = N
    else
        riemann_boundary = i
    end if
    return
end function riemann_boundary


subroutine SpaceDis(Nx, Ny, input_rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, output_rhs)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(dim*k*k,dim*k*k) :: BigMat1, BigMat2, BigMat3, BigMat4, C=0.0, F=0.0
    real :: input_rho(dim*k*k,Nx*Ny), output_rhs(dim*k*k,Nx*Ny)
    real :: x(Nx+1), y(Ny+1), Gauss(k), wt(k)
    real, allocatable :: f1(:,:), f2(:,:), FFlux1(:,:), FFlux2(:,:)
    real, allocatable :: q1(:,:), q2(:,:), QFlux1(:,:), QFlux2(:,:), RhoFlux1(:,:), RhoFlux2(:,:)
    real, allocatable :: xi(:,:), XiFlux1(:,:), XiFlux2(:,:)
    real, allocatable :: theta1(:,:), theta2(:,:)
    real :: para1, para2, g=0.001
    integer :: i,j
    
    para1 = 2.0D0*Nx/(xr-xl)
    para2 = 2.0D0*Ny/(yu-yb)

    ! problem 1
    !allocate(q1(dim*k*k,Nx*Ny), q2(dim*k*k,Nx*Ny), QFlux1(dim*k*k,Nx*Ny), QFlux2(dim*k*k,Nx*Ny))   
    !call func_G(Nx, Ny, input_rho, q1)
    !call func_G(Nx, Ny, input_rho, q2)
    !q1 = -q1
    !q2 = -q2
    !
    !call flux_AL_f(Nx, Ny, QFlux1, QFlux2, input_rho, q1, q2)
    !
    !output_rhs = -para1*matmul(BigMat1,q1) - para2*matmul(BigMat2,q2) 
    !deallocate(q1,q2)
    !output_rhs = output_rhs + para1*matmul(BigMat3,QFlux1) + para2*matmul(BigMat4,QFlux2)
    !deallocate(QFlux1,QFlux2)
    
    ! problem 2
    !allocate(theta1(dim*k*k,Nx*Ny), theta2(dim*k*k,Nx*Ny),RhoFlux1(dim*k*k,Nx*Ny), RhoFlux2(dim*k*k,Nx*Ny))
    !allocate(q1(dim*k*k,Nx*Ny), q2(dim*k*k,Nx*Ny), QFlux1(dim*k*k,Nx*Ny), QFlux2(dim*k*k,Nx*Ny))
    !
    !call flux_C_f(Nx, Ny, RhoFlux1, RhoFlux2, input_rho, input_rho, input_rho)
    !theta1 = -para1*matmul(BigMat1, input_rho) + para1*matmul(BigMat3, RhoFlux1)
    !theta2 = -para2*matmul(BigMat2, input_rho) + para2*matmul(BigMat4, RhoFlux2)
    !
    !call func_Q(Nx, Ny, theta1, input_rho, q1)
    !call func_Q(Nx, Ny, theta2, input_rho, q2)
    !deallocate(theta1, theta2, RhoFlux1, RhoFlux2)
    !
    !call flux_LF_f(Nx, Ny, QFlux1, QFlux2, input_rho, q1, q2)
    !
    !output_rhs = -para1*matmul(BigMat1,q1) - para2*matmul(BigMat2,q2) 
    !output_rhs = output_rhs + para1*matmul(BigMat3,QFlux1) + para2*matmul(BigMat4,QFlux2)
    !deallocate(q1, q2, QFlux1,QFlux2)
    
    
    ! 带源项非线性方程精确解
    allocate(theta1(dim*k*k,Nx*Ny), theta2(dim*k*k,Nx*Ny),RhoFlux1(dim*k*k,Nx*Ny), RhoFlux2(dim*k*k,Nx*Ny))
    allocate(q1(dim*k*k,Nx*Ny), q2(dim*k*k,Nx*Ny), QFlux1(dim*k*k,Nx*Ny), QFlux2(dim*k*k,Nx*Ny))
    
    call flux_C_f(Nx, Ny, RhoFlux1, RhoFlux2, input_rho, input_rho, input_rho)
    theta1 = -para1*matmul(BigMat1, input_rho) + para1*matmul(BigMat3, RhoFlux1)
    theta2 = -para2*matmul(BigMat2, input_rho) + para2*matmul(BigMat4, RhoFlux2)
    
    call func_Q(Nx, Ny, theta1, input_rho, q1)
    call func_Q(Nx, Ny, theta2, input_rho, q2)
    deallocate(theta1, theta2, RhoFlux1, RhoFlux2)
    
    !call flux_C_f(Nx, Ny, QFlux1, QFlux2, input_rho, q1, q2)
    call flux_AL_f(Nx, Ny, QFlux1, QFlux2, input_rho, q1, q2)
    
    output_rhs = -para1*matmul(BigMat1,q1) - para2*matmul(BigMat2,q2) 
    output_rhs = output_rhs + para1*matmul(BigMat3,QFlux1) + para2*matmul(BigMat4,QFlux2)
    deallocate(q1,q2,QFlux1,QFlux2)
    
    
    return
end subroutine SpaceDis
    
subroutine matrixC(rho, matC)
    use global
    implicit none
    real :: rho(dim*k*k), matC(dim*k*k,dim*k*k), rho1, rho2
    integer :: i
    do i = 1,k*k
        rho1 = rho(dim*i-1)
        rho2 = rho(dim*i)
        matC(dim*i-1,dim*i-1) = rho2 + 2.0*rho1
        matC(dim*i-1,dim*i) = rho2
        matC(dim*i,dim*i-1) = rho1
        matC(dim*i,dim*i) = rho1 +2.0*rho2
    enddo
    return
end subroutine matrixC
    
subroutine matrixF(rho, matF)
    use global
    implicit none
    real :: rho(dim*k*k), matF(dim*k*k,dim*k*k), rho1, rho2
    integer :: i
    !matF = 0.0
    do i = 1,k*k
        rho1 = rho(dim*i-1)
        rho2 = rho(dim*i)
        !SKT population problem
        !matF(dim*i-1,dim*i-1) = rho1*rho2 + 2.0*rho1**2
        !matF(dim*i-1,dim*i) = rho1*rho2
        !matF(dim*i,dim*i-1) = rho1*rho2
        !matF(dim*i,dim*i) = rho1*rho2 +2.0*rho2**2
        !Surfactant problem
        matF(dim*i-1,dim*i-1) = rho1**3/3.0D0
        matF(dim*i-1,dim*i) = 0.5D0*rho1*rho1*rho2
        matF(dim*i,dim*i-1) = 0.5D0*rho1*rho1*rho2
        matF(dim*i,dim*i) = rho1*rho2**2
    enddo
    return
end subroutine matrixF    
    
    
subroutine Error(Nx, Ny, rho, rho_exact, wt, err1, err2)
    use global
    implicit none
    integer :: Nx, Ny, x1(k*k), x2(k*k)
    real :: rho(dim*k*k,Nx*Ny), rho_exact(dim*k*k,Nx*Ny), wt(k*k)
    real :: err1(3), err2(3)
    integer :: i, j
    
    err1=0.0
    err2=0.0
    do i=1,k*k
      x1(i) = 2*(i-1) + 1
      x2(i) = 2*(i-1) + 2
    enddo
    do i = 1,Nx*Ny
        err1(1) = err1(1) + dot_product(wt, abs(rho(x1,i) - rho_exact(x1,i)))
        err2(1) = err2(1) + dot_product(wt, abs(rho(x2,i) - rho_exact(x2,i)))
        err1(2) = err1(2) + dot_product(wt, (rho(x1,i)-rho_exact(x1,i))**2)
        err2(2) = err2(2) + dot_product(wt, (rho(x2,i)-rho_exact(x2,i))**2)
    end do
    err1(1) = err1(1)/Nx/Ny
    err2(1) = err2(1)/Nx/Ny
    err1(2) = sqrt(err1(2)/Nx/Ny)
    err2(2) = sqrt(err2(2)/Nx/Ny)
    err1(3) = maxval(abs(rho(x1,:)-rho_exact(x1,:)))
    err2(3) = maxval(abs(rho(x2,:)-rho_exact(x2,:)))
    return
end subroutine Error
    
program main
    use global
    
    implicit none
   
    ! 设置参数
    integer :: i, j, Nx, Ny, LoopIter
    integer :: iter = 0
    real :: hx, hy, dx, dy, time = 0, dt, tau
    real :: local_x(k), local_y(k)
    
    real, allocatable :: Gauss(:), wt(:), longwt(:), x(:), y(:)
    real, allocatable :: rho(:,:), rho_exact(:,:), rhs(:,:)
    real, dimension(dim,dim) ::DimE=0
    real, dimension(k,k) :: D=0, E=0
    real, dimension(k*k,k*k) :: D1=0, D2=0, B1=0, B2=0, M=0, invM=0
    real, dimension(k*k,k*k) :: mat1=0.0, mat2=0.0, mat3=0.0, mat4=0.0
    real, dimension(k*k*dim, dim*k*k) :: BigMat1=0.0, BigMat2=0.0, BigMat3=0.0, BigMat4=0.0
    real, dimension(Loop, 3*2) :: errortable1=0.0D0, errortable2=0.0D0
    real, dimension(3) :: err1=0.0D0, err2=0.0D0
    real time_begin , time_end
    character(len = 30) :: filename1 = "result1.dat"
    character(len = 30) :: filename2 = "result2.dat"
    character(len = 30) :: filename_ex = "result_exact.dat"
    character(len = 30) :: filename_err = "result_error.dat"
    integer, parameter :: fileid1 = 11
    integer, parameter :: fileid2 = 12
    integer, parameter :: fileid_ex = 9
    integer, parameter :: fileid_err = 10
    character( len = 3 ) :: cTemp1, cTemp2
    
    ! 根据Gauss点数目给出Gauss点和权重
    allocate(Gauss(k), wt(k),longwt(k*k))
    call GaussLobatto_1d(k,Gauss,wt) ! x y方向公用一套Gauss积分点及权重
    
    ! 构造需要用到的矩阵
    call PartialLagrangian(Gauss, k, D)
    do i = 1,k
        E(i,i) = 1.0D0
    end do
    do i = 1,dim
        DimE(i,i) = 1.0D0
    end do
    call Kronecker(k,k,D,k,k,E,D2)
    call Kronecker(k,k,E,k,k,D,D1)
    do j = 1,k
    do i = 1,k
        M((j-1)*k+i,(j-1)*k+i) = wt(j)*wt(i)
        invM((j-1)*k+i,(j-1)*k+i) = 1.0/(wt(j)*wt(i))
        longwt((j-1)*k+i) = wt(i)*wt(j)
        if (j.EQ.1 .or. j.eq.k) then
            mat4((j-1)*k+i,(j-1)*k+i) = wt(i)
        else
            mat4((j-1)*k+i,(j-1)*k+i) = 0.0D0
        end if
        
        if (i.EQ.1 .or. i.eq.k) then
            mat3((j-1)*k+i,(j-1)*k+i) = wt(j)
        else
            mat3((j-1)*k+i,(j-1)*k+i) = 0.0D0
        end if
    end do
    end do
    mat1 = matmul(invM,matmul(transpose(D1),M))
    mat2 = matmul(invM,matmul(transpose(D2),M))
    mat3 = matmul(invM,mat3)
    mat4 = matmul(invM,mat4)
    
    call Kronecker(k*k,k*k,mat1,dim,dim,DimE,BigMat1)
    call Kronecker(k*k,k*k,mat2,dim,dim,DimE,BigMat2)
    call Kronecker(k*k,k*k,mat3,dim,dim,DimE,BigMat3)
    call Kronecker(k*k,k*k,mat4,dim,dim,DimE,BigMat4)
    
    call CreateFluxPara
    call CreateDomainPara
    
    ! 循环
    do LoopIter = 1,Loop
    
    ! 网格点坐标
    !Nx = 20*2**(LoopIter-1)    ! x方向小区间数量
    !Ny = 20*2**(LoopIter-1)    ! y方向小区间数量
    Nx = 20*LoopIter
    Ny = 20*LoopIter
    hx = (xr-xl)/Nx
    hy = (yu-yb)/Ny
    tau = 0.00025D0*((hx+hy)/2.0D0)**2     ! 固定时间步长问题
    
    allocate(x(Nx+1), y(Ny+1))
    
    do i = 1,Nx+1
        x(i) = xl + hx*(i-1)    
    end do
    do j = 1,Ny+1
        y(j) = yb + hy*(j-1)
    end do
    
    ! 设定初值
    allocate(rho(dim*k*k,Nx*Ny),rho_exact(dim*k*k,Nx*Ny),rhs(dim*k*k,Nx*Ny))
    call initial(Nx, Ny, rho, Gauss, x, y)
    
    ! 主循环
    call CPU_TIME(time_begin)
    rhs = 0.D0
    time = 0.0D0
    do while(time+tau < T)
        call SSP_RK(Nx, Ny, rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs, tau, time)
        time = time + tau
        iter = iter + 1
        !write(*,*) "time=", time
    end do
    dt = T - time
    call SSP_RK(Nx, Ny, rho, x, y, Gauss, wt, BigMat1, BigMat2, BigMat3, BigMat4, rhs, dt, time)
    time = time + dt
    call CPU_TIME(time_end)

    ! 精确解
    call exact(Nx, Ny, rho_exact, Gauss, x, y, time)

    ! 计算误差
    call Error(Nx, Ny, rho, rho_exact, longwt, err1, err2)
    write(*,*) "Grid size", Nx, Ny
    write(*,*) "error norm"
    write(*,'(3(E12.4))') err1(1), err1(2), err1(3)
    write(*,'(3(E12.4))') err2(1), err2(2), err2(3)
    
    write(*,*) "iteration", iter
    write(*,*) "ending time", time
    write(*,*) "time cost", time_end - time_begin
    
    errortable1(LoopIter,1) = err1(1)
    errortable1(LoopIter,3) = err1(2)
    errortable1(LoopIter,5) = err1(3)
    errortable2(LoopIter,1) = err2(1)
    errortable2(LoopIter,3) = err2(2)
    errortable2(LoopIter,5) = err2(3)
    
    ! 将数据输出到txt文件
    !open(unit = fileid1, file = filename1)
    !open(unit = fileid2, file = filename2)
    !write(fileid1,*) 'VARIABLES = "X","Y","Z" '
    !write(fileid1,*) 'zone ', 'I= ',Nx, 'J= ',Ny, ' F=point' 
    !write(fileid2,*) 'VARIABLES = "X","Y","Z" '
    !write(fileid2,*) 'zone ', 'I= ',Nx, 'J= ',Ny, ' F=point' 
    !do j=1,Nx
    !do i=1,Ny
    !    write(fileid1,*) (x(i)+x(i+1))/2, (y(j)+y(j+1))/2, rho(9,(j-1)*Nx+i)
    !    write(fileid2,*) (x(i)+x(i+1))/2, (y(j)+y(j+1))/2, rho(10,(j-1)*Nx+i)
    !enddo
    !enddo
    !close(fileid1)
    !close(fileid2)
    
    
    !write(cTemp1,'(i3)') Nx
    !write(cTemp2,'(i3)') Ny
    !filename = 'Acc_Nx' // trim(adjustl(cTemp1)) // 'Ny' // trim(adjustl(cTemp2)) // '_rho.dat'
    !open(unit = fileid1, file = filename1)
    !open(unit = fileid2, file = filename2)
    !open(unit = fileid_ex, file = filename_ex)
    !open(unit = fileid_err, file = filename_err)
    !write(fileid1,*) 'VARIABLES = "X","Y","Z" '
    !write(fileid1,*) 'zone ', 'I= ',Nx, 'J= ',Ny, ' F=point' 
    !write(fileid2,*) 'VARIABLES = "X","Y","Z" '
    !write(fileid2,*) 'zone ', 'I= ',Nx, 'J= ',Ny, ' F=point' 
    !write(fileid_ex,*) 'VARIABLES = "X","Y","Z" '
    !write(fileid_ex,*) 'zone ', 'I= ',Nx, 'J= ',Ny, ' F=point' 
    !write(fileid_err,*) 'VARIABLES = "X","Y","Z" '
    !write(fileid_err,*) 'zone ', 'I= ',Nx, 'J= ',Ny, ' F=point' 
    !do j=1,Nx
    !do i=1,Ny
    !    write(fileid1,*) (x(i)+x(i+1))/2, (y(j)+y(j+1))/2, rho(9,(j-1)*Nx+i)
    !    write(fileid2,*) (x(i)+x(i+1))/2, (y(j)+y(j+1))/2, rho(10,(j-1)*Nx+i)
    !    write(fileid_ex,*) (x(i)+x(i+1))/2, (y(j)+y(j+1))/2, rho_exact(10,(j-1)*Nx+i)
    !    if (abs(rho(5,(j-1)*Nx+i)-rho_exact(5,(j-1)*Nx+i)) < 1e-8) then
    !        write(fileid_err,*) (x(i)+x(i+1))/2, (y(j)+y(j+1))/2, -16.0D0
    !    else
    !        write(fileid_err,*) (x(i)+x(i+1))/2, (y(j)+y(j+1))/2, abs(rho(10,(j-1)*Nx+i)-rho_exact(10,(j-1)*Nx+i))
    !    end if
    !!    100 format(3(f20.16, ','), ';')
    !enddo
    !enddo
    !close(fileid1)
    !close(fileid2)
    !close(fileid_ex)
    !close(fileid_err)

    deallocate(x,y,rho,rho_exact,rhs)

    enddo   ! 结束循环Loop
    
    do i = 1,Loop-1
        !errortable1(i+1,2) = LOG(errortable1(i,1)/errortable1(i+1,1))/LOG(2.0D0)
        !errortable1(i+1,4) = LOG(errortable1(i,3)/errortable1(i+1,3))/LOG(2.0D0)
        !errortable1(i+1,6) = LOG(errortable1(i,5)/errortable1(i+1,5))/LOG(2.0D0)
        !
        !errortable2(i+1,2) = LOG(errortable2(i,1)/errortable2(i+1,1))/LOG(2.0D0)
        !errortable2(i+1,4) = LOG(errortable2(i,3)/errortable2(i+1,3))/LOG(2.0D0)
        !errortable2(i+1,6) = LOG(errortable2(i,5)/errortable2(i+1,5))/LOG(2.0D0)
        
        errortable1(i+1,2) = LOG(errortable1(i,1)/errortable1(i+1,1))/LOG(1.0D0*(i+1)/i)
        errortable1(i+1,4) = LOG(errortable1(i,3)/errortable1(i+1,3))/LOG(1.0D0*(i+1)/i)
        errortable1(i+1,6) = LOG(errortable1(i,5)/errortable1(i+1,5))/LOG(1.0D0*(i+1)/i)
        
        errortable2(i+1,2) = LOG(errortable2(i,1)/errortable2(i+1,1))/LOG(1.0D0*(i+1)/i)
        errortable2(i+1,4) = LOG(errortable2(i,3)/errortable2(i+1,3))/LOG(1.0D0*(i+1)/i)
        errortable2(i+1,6) = LOG(errortable2(i,5)/errortable2(i+1,5))/LOG(1.0D0*(i+1)/i)
    enddo
    
    write(*,*) "Rho1 Error Norm"
    do i = 1,Loop
        write(*,200) errortable1(i,:)        
    enddo
    
    write(*,*) "Rho2 Error Norm"
    do i = 1,Loop
        write(*,200) errortable2(i,:)     
    enddo
    200 format(E11.3, F6.2, E11.3, F6.2, E11.3, F6.2)   
    
    deallocate(Gauss,wt,longwt)
    call DeleteFluxPara

      
stop
end program main