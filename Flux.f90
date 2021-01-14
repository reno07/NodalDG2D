subroutine Flux_DW_f(Nx, Ny, FFluxX, FFluxY, rho, f1, f2)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(dim*k*k, Nx*Ny) :: FFluxX, FFluxY, rho, f1, f2
    integer :: i, j, ii, jj, l, kk, fkk, Line, XLineFlux, YLineFlux
    integer, external :: periodic_boundary
    real :: alpha
    
    ! 对于2D下每一个单元 x y方向数值通量交互 单元左右两条边界 
    ! 这些点在该单元索引为 BoundaryIndex
    ! 对应在相邻单元为 ii,jj 行号为LineFlux 点的索引为XFluxIndex
        
    ! Downwind flux
    do j = 1,Ny
        do i = 1,Nx
            Line = (j-1)*Nx + i
            do l = 1,4*k-4
            ii = periodic_boundary(i+n1(l),Nx) 
            jj = periodic_boundary(j+n2(l),Ny)
            kk = BoundaryIndex(l)
            XLineFlux = (j-1)*Nx + ii
            YLineFlux = (jj-1)*Nx + i
            
            if (n1(l) < 0) then
                FFluxX(kk*dim-1,Line) = -f1(kk*dim-1,Line)
                FFluxX(kk*dim,Line) = -f1(kk*dim,Line)
            else if (n1(l) > 0) then
                FFluxX(kk*dim-1,Line) = f1(dim*XFluxIndex(l)-1,XLineFlux)
                FFluxX(kk*dim,Line) = f1(dim*XFluxIndex(l),XLineFlux)
            else
                FFluxX(kk*dim-1,Line) = 0.0D0
                FFluxX(kk*dim,Line) = 0.0D0
            endif
            
            if (n2(l) < 0) then
                FFluxY(dim*kk-1,Line) = -f2(dim*kk-1,Line)
                FFluxY(dim*kk,Line) = -f2(dim*kk,Line)
            else if (n2(l) > 0) then
                FFluxY(kk*dim-1,Line) = f2(dim*YFluxIndex(l)-1,YLineFlux)
                FFluxY(kk*dim,Line) = f2(dim*YFluxIndex(l),YLineFlux)
            else
                FFluxY(dim*kk-1,Line) = 0.0D0
                FFluxY(dim*kk,Line) = 0.0D0
            endif
            
            enddo
        enddo
    enddo
    return
end subroutine flux_DW_f  
    
subroutine Flux_UW_f(Nx, Ny, FFluxX, FFluxY, rho, f1, f2)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(dim*k*k, Nx*Ny) :: FFluxX, FFluxY, rho, f1, f2
    integer :: i, j, ii, jj, l, kk, fkk, Line, XLineFlux, YLineFlux
    integer, external :: periodic_boundary, riemann_boundary
    real :: alpha
    
    ! 对于2D下每一个单元 x方向数值通量交互 单元左右两条边界 
    ! 这些点在该单元索引为 BoundaryIndex
    ! 对应在相邻单元为 ii,jj 行号为LineFlux 点的索引为XFluxIndex

    ! Upwind flux
    do j = 1,Ny
        do i = 1,Nx
            Line = (j-1)*Nx + i
            do l = 1,4*k-4
            ii = periodic_boundary(i+n1(l),Nx) 
            jj = periodic_boundary(j+n2(l),Ny)
            
            kk = BoundaryIndex(l)
            XLineFlux = (j-1)*Nx + ii
            YLineFlux = (jj-1)*Nx + i
            
            if (n1(l) > 0) then
                FFluxX(kk*dim-1,Line) = f1(kk*dim-1,Line)
                FFluxX(kk*dim,Line) = f1(kk*dim,Line)
            else if (n1(l) < 0) then
                FFluxX(kk*dim-1,Line) = -f1(dim*XFluxIndex(l)-1,XLineFlux)
                FFluxX(kk*dim,Line) = -f1(dim*XFluxIndex(l),XLineFlux)
            else
                FFluxX(kk*dim-1,Line) = 0.0D0
                FFluxX(kk*dim,Line) = 0.0D0
            endif
            
            if (n2(l) > 0) then
                FFluxY(dim*kk-1,Line) = f2(dim*kk-1,Line)
                FFluxY(dim*kk,Line) = f2(dim*kk,Line)
            else if (n2(l) < 0) then
                FFluxY(kk*dim-1,Line) = -f2(dim*YFluxIndex(l)-1,YLineFlux)
                FFluxY(kk*dim,Line) = -f2(dim*YFluxIndex(l),YLineFlux)
            else
                FFluxY(dim*kk-1,Line) = 0.0D0
                FFluxY(dim*kk,Line) = 0.0D0
            endif
            
            enddo
        enddo
    enddo
    return
end subroutine flux_UW_f  
    
subroutine Flux_LF_f(Nx, Ny, FFluxX, FFluxY, rho, f1, f2)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(dim*k*k, Nx*Ny) :: FFluxX, FFluxY, rho, f1, f2
    integer :: i, j, ii, jj, l, kk, fkk, Line, XLineFlux, YLineFlux
    integer, external :: periodic_boundary
    real :: alpha1, alpha2
    
    ! 对于2D下每一个单元 x方向数值通量交互 单元左右两条边界 
    ! 这些点在该单元索引为 BoundaryIndex
    ! 对应在相邻单元为 ii,jj 行号为LineFlux 点的索引为XFluxIndex
        
    ! Lax-Friedrichs flux
    do j = 1,Ny
        do i = 1,Nx
            Line = (j-1)*Nx + i
            do l = 1,4*k-4
            ii = periodic_boundary(i+n1(l),Nx) 
            jj = periodic_boundary(j+n2(l),Ny)
            kk = BoundaryIndex(l)
            XLineFlux = (j-1)*Nx + ii
            YLineFlux = (jj-1)*Nx + i
            
            alpha1 = max(abs(f1(kk*dim-1,Line)/rho(kk*dim-1,Line)), abs(f1(XFluxIndex(l)*dim-1,XLineFlux)/rho(XFluxIndex(l)*dim-1,XLineFlux)))/8.0
            alpha2 = max(abs(f1(kk*dim,Line)/rho(kk*dim,Line)), abs(f1(XFluxIndex(l)*dim,XLineFlux)/rho(XFluxIndex(l)*dim,XLineFlux)))/8.0   
            
            !alpha1 = max(rho, 1.0)
            if (n1(l) > 0) then
                FFluxX(kk*dim-1,Line) = 0.5D0*(f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux)) + alpha1*(rho(XFluxIndex(l)*dim-1,XLineFlux) - rho(kk*dim-1,Line))
                FFluxX(kk*dim,Line) = 0.5D0*(f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux)) + alpha2*(rho(XFluxIndex(l)*dim,XLineFlux) - rho(kk*dim,Line))
            else if (n1(l) < 0) then
                FFluxX(kk*dim-1,Line) = -0.5D0*(f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux)) - alpha1*(rho(kk*dim-1,Line) - rho(XFluxIndex(l)*dim-1,XLineFlux))
                FFluxX(kk*dim,Line) = -0.5D0*(f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux)) - alpha2*(rho(kk*dim,Line) - rho(XFluxIndex(l)*dim,XLineFlux))
            else
                FFluxX(kk*dim-1,Line) = 0.0D0
                FFluxX(kk*dim,Line) = 0.0D0
            endif
            
            alpha1 = max(abs(f2(kk*dim-1,Line)/rho(kk*dim-1,Line)), abs(f2(YFluxIndex(l)*dim-1,YLineFlux)/rho(YFluxIndex(l)*dim-1,YLineFlux)))/8.0
            alpha2 = max(abs(f2(kk*dim,Line)/rho(kk*dim,Line)), abs(f2(YFluxIndex(l)*dim,YLineFlux)/rho(YFluxIndex(l)*dim,YLineFlux)))/8.0        
            if (n2(l) > 0) then
                FFluxY(kk*dim-1,Line) = 0.5D0*(f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux)) + alpha1*(rho(YFluxIndex(l)*dim-1,YLineFlux) - rho(kk*dim-1,Line))
                FFluxY(kk*dim,Line) = 0.5D0*(f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux)) + alpha2*(rho(YFluxIndex(l)*dim,YLineFlux) - rho(kk*dim,Line))
            else if (n2(l) < 0) then
                FFluxY(kk*dim-1,Line) = -0.5D0*(f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux)) - alpha1*(rho(kk*dim-1,Line) - rho(YFluxIndex(l)*dim-1,YLineFlux))
                FFluxY(kk*dim,Line) = -0.5D0*(f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux)) - alpha2*(rho(kk*dim,Line) - rho(YFluxIndex(l)*dim,YLineFlux))
            else
                FFluxY(kk*dim-1,Line) = 0.0D0
                FFluxY(kk*dim,Line) = 0.0D0
            endif
            
            enddo
        enddo
    enddo
    
    return
end subroutine flux_LF_f   
    
subroutine Flux_C_f(Nx, Ny, FFluxX, FFluxY, rho, f1, f2)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(dim*k*k, Nx*Ny) :: FFluxX, FFluxY, rho, f1, f2
    integer :: i, j, ii, jj, l, kk, fkk, Line, XLineFlux, YLineFlux
    integer, external :: periodic_boundary
    real :: alpha1, alpha2
    
    ! 对于2D下每一个单元 x方向数值通量交互 单元左右两条边界 
    ! 这些点在该单元索引为 BoundaryIndex
    ! 对应在相邻单元为 ii,jj 行号为LineFlux 点的索引为XFluxIndex
        
    ! Central flux
    do j = 1,Ny
        do i = 1,Nx
            Line = (j-1)*Nx + i
            do l = 1,4*k-4
            ii = periodic_boundary(i+n1(l),Nx) 
            jj = periodic_boundary(j+n2(l),Ny)
            kk = BoundaryIndex(l)
            XLineFlux = (j-1)*Nx + ii
            YLineFlux = (jj-1)*Nx + i

            if (n1(l) > 0) then
                FFluxX(kk*dim-1,Line) = 0.5D0*(f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux)) 
                FFluxX(kk*dim,Line) = 0.5D0*(f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux))
            else if (n1(l) < 0) then
                FFluxX(kk*dim-1,Line) = -0.5D0*(f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux)) 
                FFluxX(kk*dim,Line) = -0.5D0*(f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux))
            else
                FFluxX(kk*dim-1,Line) = 0.0D0
                FFluxX(kk*dim,Line) = 0.0D0
            endif
            
            if (n2(l) > 0) then
                FFluxY(kk*dim-1,Line) = 0.5D0*(f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux))
                FFluxY(kk*dim,Line) = 0.5D0*(f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux)) 
            else if (n2(l) < 0) then
                FFluxY(kk*dim-1,Line) = -0.5D0*(f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux))
                FFluxY(kk*dim,Line) = -0.5D0*(f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux))
            else
                FFluxY(kk*dim-1,Line) = 0.0D0
                FFluxY(kk*dim,Line) = 0.0D0
            endif
            
            enddo
        enddo
    enddo
    return
end subroutine flux_C_f     
    
subroutine Flux_ZeroUW_f(Nx, Ny, FFluxX, FFluxY, rho, f1, f2)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(dim*k*k, Nx*Ny) :: FFluxX, FFluxY, rho, f1, f2
    integer :: i, j, ii, jj, l, kk, fkk, Line, XLineFlux, YLineFlux
    integer, external :: periodic_boundary, riemann_boundary
    real :: alpha
    
    ! zero-flux 边界条件 中间取upwind
    do j = 1,Ny
        do i = 1,Nx
            Line = (j-1)*Nx + i
            do l = 1,4*k-4
            kk = BoundaryIndex(l)
            if (i.EQ.1 .or. i.EQ.Nx .or. j.EQ.1 .or. j.EQ.Ny) then
                FFluxX(kk*dim-1,Line) = 0.0D0
                FFluxX(kk*dim,Line) = 0.0D0
                FFluxY(dim*kk-1,Line) = 0.0D0
                FFluxY(dim*kk,Line) = 0.0D0
            else
                ii = riemann_boundary(i+n1(l),Nx) 
                jj = riemann_boundary(j+n2(l),Ny)
                XLineFlux = (j-1)*Nx + ii
                YLineFlux = (jj-1)*Nx + i
            
                if (n1(l) > 0) then
                    FFluxX(kk*dim-1,Line) = f1(kk*dim-1,Line)
                    FFluxX(kk*dim,Line) = f1(kk*dim,Line)
                else if (n1(l) < 0) then
                    FFluxX(kk*dim-1,Line) = -f1(dim*XFluxIndex(l)-1,XLineFlux)
                    FFluxX(kk*dim,Line) = -f1(dim*XFluxIndex(l),XLineFlux)
                else
                    FFluxX(kk*dim-1,Line) = 0.0D0
                    FFluxX(kk*dim,Line) = 0.0D0
                endif
            
                if (n2(l) > 0) then
                    FFluxY(dim*kk-1,Line) = f2(dim*kk-1,Line)
                    FFluxY(dim*kk,Line) = f2(dim*kk,Line)
                else if (n2(l) < 0) then
                    FFluxY(kk*dim-1,Line) = -f2(dim*YFluxIndex(l)-1,YLineFlux)
                    FFluxY(kk*dim,Line) = -f2(dim*YFluxIndex(l),YLineFlux)
                else
                    FFluxY(dim*kk-1,Line) = 0.0D0
                    FFluxY(dim*kk,Line) = 0.0D0
                endif
            endif
            
            enddo
        enddo
    enddo
    return
end subroutine flux_ZeroUW_f      
    
    
subroutine Flux_ZeroC_f(Nx, Ny, FFluxX, FFluxY, rho, f1, f2)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(dim*k*k, Nx*Ny) :: FFluxX, FFluxY, rho, f1, f2
    integer :: i, j, ii, jj, l, kk, fkk, Line, XLineFlux, YLineFlux
    integer, external :: periodic_boundary
    real :: alpha1, alpha2
    
    ! 对于2D下每一个单元 x方向数值通量交互 单元左右两条边界 
    ! 这些点在该单元索引为 BoundaryIndex
    ! 对应在相邻单元为 ii,jj 行号为LineFlux 点的索引为XFluxIndex
        
    ! Central flux
    do j = 1,Ny
        do i = 1,Nx
            Line = (j-1)*Nx + i
            do l = 1,4*k-4
            kk = BoundaryIndex(l)
            if (i.EQ.1 .or. i.EQ.Nx) then
                FFluxX(kk*dim-1,Line) = 0.0D0
                FFluxX(kk*dim,Line) = 0.0D0
            else if (j.EQ.1 .or. j.EQ.Ny) then
                FFluxY(dim*kk-1,Line) = 0.0D0
                FFluxY(dim*kk,Line) = 0.0D0
            else
                ii = periodic_boundary(i+n1(l),Nx) 
                jj = periodic_boundary(j+n2(l),Ny)
                XLineFlux = (j-1)*Nx + ii
                YLineFlux = (jj-1)*Nx + i
            
                alpha1 = max(abs(rho(kk*dim-1,Line)), abs(rho(XFluxIndex(l)*dim-1,XLineFlux)))/4.0
                alpha2 = max(abs(rho(kk*dim,Line)), abs(rho(XFluxIndex(l)*dim,XLineFlux)))/4.0        
                alpha1 = 0.0D0
                alpha2 = 0.0D0
                if (n1(l) > 0) then
                    FFluxX(kk*dim-1,Line) = 0.5D0*(f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux)) - alpha1*(rho(XFluxIndex(l)*dim-1,XLineFlux) - rho(kk*dim-1,Line))
                    FFluxX(kk*dim,Line) = 0.5D0*(f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux)) - alpha1*(rho(XFluxIndex(l)*dim,XLineFlux) - rho(kk*dim,Line))
                else if (n1(l) < 0) then
                    FFluxX(kk*dim-1,Line) = -0.5D0*(f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux)) + alpha1*(rho(kk*dim-1,Line) - rho(XFluxIndex(l)*dim-1,XLineFlux))
                    FFluxX(kk*dim,Line) = -0.5D0*(f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux)) + alpha1*(rho(kk*dim,Line) - rho(XFluxIndex(l)*dim,XLineFlux))
                else
                    FFluxX(kk*dim-1,Line) = 0.0D0
                    FFluxX(kk*dim,Line) = 0.0D0
                endif
            
                alpha1 = max(abs(rho(kk*dim-1,Line)), abs(rho(YFluxIndex(l)*dim-1,YLineFlux)))/4.0
                alpha2 = max(abs(rho(kk*dim,Line)), abs(rho(YFluxIndex(l)*dim,YLineFlux)))/4.0        
                alpha1 = 0.0D0
                alpha2 = 0.0D0
                if (n2(l) > 0) then
                    FFluxY(kk*dim-1,Line) = 0.5D0*(f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux)) - alpha2*(rho(YFluxIndex(l)*dim-1,YLineFlux) - rho(kk*dim-1,Line))
                    FFluxY(kk*dim,Line) = 0.5D0*(f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux)) - alpha2*(rho(YFluxIndex(l)*dim,YLineFlux) - rho(kk*dim,Line))
                else if (n2(l) < 0) then
                    FFluxY(kk*dim-1,Line) = -0.5D0*(f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux)) + alpha2*(rho(kk*dim-1,Line) - rho(YFluxIndex(l)*dim-1,YLineFlux))
                    FFluxY(kk*dim,Line) = -0.5D0*(f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux)) + alpha2*(rho(kk*dim,Line) - rho(YFluxIndex(l)*dim,YLineFlux))
                else
                    FFluxY(kk*dim-1,Line) = 0.0D0
                    FFluxY(kk*dim,Line) = 0.0D0
                endif
            endif
            
            enddo
        enddo
    enddo
    return
end subroutine flux_ZeroC_f        
    
    
subroutine Flux_AL_f(Nx, Ny, FFluxX, FFluxY, rho, f1, f2)
    use global
    implicit none
    integer :: Nx, Ny
    real, dimension(dim*k*k, Nx*Ny) :: FFluxX, FFluxY, rho, f1, f2
    integer :: i, j, ii, jj, l, kk, fkk, Line, XLineFlux, YLineFlux
    integer, external :: periodic_boundary
    real :: alpha1, alpha2
    
    ! 对于2D下每一个单元 x方向数值通量交互 单元左右两条边界 
    ! 这些点在该单元索引为 BoundaryIndex
    ! 对应在相邻单元为 ii,jj 行号为LineFlux 点的索引为XFluxIndex
        
    ! Alternative flux
    do j = 1,Ny
        do i = 1,Nx
            Line = (j-1)*Nx + i
            do l = 1,4*k-4
            ii = periodic_boundary(i+n1(l),Nx) 
            jj = periodic_boundary(j+n2(l),Ny)
            kk = BoundaryIndex(l)
            XLineFlux = (j-1)*Nx + ii
            YLineFlux = (jj-1)*Nx + i
                  
            alpha1 = -0.5D0/pi
            alpha2 = -0.5D0*pi
            
            alpha1 = -0.5D0/DSQRT(abs(rho(kk*dim-1,Line)**1 + rho(XFluxIndex(l)*dim-1,XLineFlux)**1)/2.0)
            alpha2 = -0.5D0*DSQRT(abs(rho(kk*dim-1,Line)**1 + rho(XFluxIndex(l)*dim-1,XLineFlux)**1)/2.0)
            
            if (n1(l) > 0) then
                FFluxX(kk*dim-1,Line) = 0.5D0*(f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux)) &
                            + alpha1*(-f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux))
                FFluxX(kk*dim,Line) = 0.5D0*(f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux)) &
                            + alpha2*(-f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux))
            else if (n1(l) < 0) then
                FFluxX(kk*dim-1,Line) = -0.5D0*(f1(kk*dim-1,Line) + f1(XFluxIndex(l)*dim-1,XLineFlux)) &
                            - alpha1*(f1(kk*dim,Line) - f1(XFluxIndex(l)*dim,XLineFlux))
                FFluxX(kk*dim,Line) = -0.5D0*(f1(kk*dim,Line) + f1(XFluxIndex(l)*dim,XLineFlux)) &
                            - alpha2*(f1(kk*dim-1,Line) - f1(XFluxIndex(l)*dim-1,XLineFlux))
            else
                FFluxX(kk*dim-1,Line) = 0.0D0
                FFluxX(kk*dim,Line) = 0.0D0
            endif
            
            alpha1 = -0.5D0/DSQRT(abs(rho(kk*dim-1,Line)**1 + rho(YFluxIndex(l)*dim-1,YLineFlux)**1)/2.0)
            alpha2 = -0.5D0*DSQRT(abs(rho(kk*dim-1,Line)**1 + rho(YFluxIndex(l)*dim-1,YLineFlux)**1)/2.0)
      
            if (n2(l) > 0) then
                FFluxY(kk*dim-1,Line) = 0.5D0*(f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux)) &
                            + alpha1*(-f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux))
                FFluxY(kk*dim,Line) = 0.5D0*(f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux)) &
                            + alpha2*(-f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux))
            else if (n2(l) < 0) then
                FFluxY(kk*dim-1,Line) = -0.5D0*(f2(kk*dim-1,Line) + f2(YFluxIndex(l)*dim-1,YLineFlux)) &
                            - alpha1*(f2(kk*dim,Line) - f2(YFluxIndex(l)*dim,YLineFlux))
                FFluxY(kk*dim,Line) = -0.5D0*(f2(kk*dim,Line) + f2(YFluxIndex(l)*dim,YLineFlux)) &
                            - alpha2*(f2(kk*dim-1,Line) - f2(YFluxIndex(l)*dim-1,YLineFlux))
            else
                FFluxY(kk*dim-1,Line) = 0.0D0
                FFluxY(kk*dim,Line) = 0.0D0
            endif
            
            enddo
        enddo
    enddo
    return
end subroutine flux_AL_f         
    