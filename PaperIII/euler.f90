! File contains the fluxes for the two-dimensional euler equations and the corresponding entropy fluxes
module euler
        implicit none
        double precision :: gamma = 1.4D0
        contains

        ! calculates pressure from conserved variables
        function p(u)
                implicit none
                double precision, dimension(4) :: u
                double precision :: p
                p = (gamma-1)*(u(4)-0.5*(u(2)**2 + u(3)**2)/u(1))
        end function p

        ! derivative of the pressure from conserved variables
        function dp(u)
                implicit none
                double precision, dimension(4) :: u
                double precision, dimension(4) :: dp
                dp(1) = (gamma-1)*0.5*(u(2)**2 + u(3)**2)*u(1)**(-2.0)
                dp(2) = -(gamma-1)*(u(2)/u(1))
                dp(3) = -(gamma-1)*(u(3)/u(1))
                dp(4) = (gamma-1)
        end function dp

        ! calculates speed of sound from conserved variables
        function a(u)
                implicit none
                double precision, dimension(4) :: u
                double precision :: a
                a = sqrt(gamma*p(u)/ u(1))
        end function a

        ! Calcualtes the 'numeric' sound speed - used to remove problems with
        ! sonic point glitches following Tang:"On the sonic point glitch"
        ! in the HLL numerical flux.
        function anum(u)
                implicit none
                double precision, dimension(4) :: u
                double precision :: anum
                anum = sqrt(2*gamma*p(u)/ u(1))
        end function anum


        ! speed of sound directly from pressure and density
        function apr(press, rho)
                implicit none
                double precision :: press, rho
                double precision :: apr
                apr = sqrt(gamma*press/rho)
        end function apr

        ! 'numerical' speed of sound directly from pressure and density
        function aprnum(press, rho)
                implicit none
                double precision :: press, rho
                double precision :: aprnum
                aprnum = sqrt(2*gamma*press/rho)
        end function aprnum

        ! Converts between physical and conserved variables
        function ConsVar(rho, vx, vy, p)
                implicit none
                double precision, dimension(4) :: ConsVar
                double precision :: rho, vx, vy, p
                ConsVar(1) = rho
                ConsVar(2) = rho*vx
                ConsVar(3) = rho*vy
                ConsVar(4) = p / (gamma-1) + 0.5 * (vx**2 + vy**2)*rho
        end function ConsVar

        ! Converts from conserved variables to physical variables
        function PhysVar(rho, rhovx, rhovy, E)
                implicit none
                double precision, dimension(4) :: PhysVar
                double precision :: rho, rhovx, rhovy, E
                PhysVar(1) = rho
                PhysVar(2) = rhovx/rho
                PhysVar(3) = rhovy/rho
                PhysVar(4) = (gamma-1)*(E - 0.5* (rhovx**2 + rhovy**2) / rho)
        end function PhysVar


        ! calculates x component of the Euler flux
        function fEuler(u)
                implicit none
                double precision, dimension(4) :: u
                double precision, dimension(4) :: fEuler
                double precision :: pu
                pu = p(u)
                fEuler(1) = u(2)
                fEuler(2) = u(2)**2/u(1) + pu
                fEuler(3) = u(2)*u(3)/u(1)
                fEuler(4) = u(2)/u(1)*(u(4) + pu)
        end function fEuler

        ! calculates y component of Euler flux
        function gEuler(u)
                implicit none
                double precision, dimension(4) :: u
                double precision, dimension(4) :: gEuler
                double precision :: pu
                pu = p(u)
                gEuler(1) = u(3)
                gEuler(2) = u(2)*u(3)/u(1)
                gEuler(3) = u(3)**2/u(1) + pu
                gEuler(4) = u(3)/u(1)*(u(4) + pu)
        end function gEuler

        ! calculates Lax-Friedrichs flux in x
        function fEulerLF(ui, uo, n, lambda)
                double precision, dimension(4) :: fEulerLF
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                double precision, intent(in) :: lambda
                fEulerLF = (fEuler(ui) + fEuler(uo) + sign(1.0D0, n(1))*(ui - uo)/(2*lambda))/2
        end function fEulerLF

        ! Lax-Friedrichs flux in y
        function gEulerLF(ui, uo, n, lambda)
                double precision, dimension(4) :: gEulerLF
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                double precision, intent(in) :: lambda
                gEulerLF = (gEuler(ui) + gEuler(uo) + sign(1.0D0, n(2))*(ui - uo)/(2*lambda))/2
        end function gEulerLF

        ! calculates the maximum speed for the conection of two given states
        function cmax(ul, ur)
                double precision, dimension(4) :: ul, ur
                double precision :: cmax
                double precision, dimension(2) :: vl, vr
                vl = ul(2:3) / ul(1)
                vr = ur(2:3) / ur(1)
                cmax = max(norm2(vl), norm2(vr)) + max(anum(ul), anum(ur))
        end function cmax


        ! Calcualtes local Lax-Friedrichs flux in x
        function fEulerLLF(ui, uo, n)
                double precision, dimension(4) :: fEulerLLF
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                fEulerLLF = (fEuler(ui) + fEuler(uo) + 1.0D0*n(1)/norm2(n)*cmax(ui, uo)*(ui - uo))/2.0D0
        end function fEulerLLF

        ! local Lax-Friedrichs flux in y
        function gEulerLLF(ui, uo, n)
                double precision, dimension(4) :: gEulerLLF
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                gEulerLLF = (gEuler(ui) + gEuler(uo) + 1.0D0*n(2)/norm2(n)*cmax(ui, uo)*(ui - uo))/2.0D0
        end function gEulerLLF

        ! Normal Lax-Friedrichs flux
        function EulerLLF(ui, uo, n)
                double precision, dimension(4) :: EulerLLF
                double precision, dimension(4) :: ui, uo
                double precision, dimension(2) :: n
                double precision, dimension(2) :: nn
                nn = n / norm2(n)
                EulerLLF = 0.5*(nn(1)*(fEuler(ui) + fEuler(uo)) + nn(2)*(gEuler(ui) + gEuler(uo)) &
                        + cmax(ui, uo)*(ui - uo))
        end function EulerLLF

        ! direct computation of the HLL flux
        function EulerHLLnf(ui, uo, n) result (EulerHLL)
                implicit none
                double precision, dimension(4) :: EulerHLL
                double precision, dimension(4) :: ui, uo
                double precision, dimension(2) :: n
                ! Auxillaries
                double precision, dimension(2) :: nn, no
                double precision, dimension(4) :: fi, fo
                double precision :: ai, ao, vi, vo, msi, mso
                nn = n / norm2(n)
                no(1) = -nn(2)
                no(2) = nn(1)

                fi = nn(1)*fEuler(ui) + nn(2)*gEuler(ui)
                fo = nn(1)*fEuler(uo) + nn(2)*gEuler(uo)
                vi = dot_product(nn, ui(2:3))/ui(1)
                vo = dot_product(nn, uo(2:3))/uo(1)
                msi = dot_product(no, ui(2:3)/ui(1))
                mso = dot_product(no, uo(2:3)/uo(1))

                ai = min(vi, vo) - max(anum(ui), anum(uo)) - abs(msi - mso)
                ao = max(vi, vo) + max(anum(ui), anum(uo)) + abs(msi - mso)

                if (0.0D0 < ai) then                                                            ! Supersonic from left
                        EulerHLL = fi
                elseif (0.0D0 < ao) then                                                        ! Transonic case
                        EulerHLL = (ao*fi - ai*fo + ai*ao*(uo - ui)) / (ao - ai)
                else                                                                            ! Supersonic from the right
                        EulerHLL = fo
                end if
        end function EulerHLLnf

        ! Computation using the rotation onto the x direction
        function EulerHLLrot(ui, uo, n) result (EulerHLL)
                implicit none
                double precision, dimension(4) :: EulerHLL
                double precision, dimension(4) :: ui, uo
                double precision, dimension(2) :: n
                ! Auxillaries
                double precision, dimension(2) :: nn, vl, vr
                double precision, dimension(4) :: fl, fr, ul, ur, um
                double precision, dimension(2, 2) :: R, Ri
                double precision :: al, ar
                nn = n / norm2(n)

                ! Rotation into the x direction
                R(1, :) = nn
                R(2, 1) = -nn(2)
                R(2, 2) = nn(1)
                ! Rotation Back
                Ri(1, 1) = nn(1)
                Ri(1, 2) = -nn(2)
                Ri(2, 1) = nn(2)
                Ri(2, 2) = nn(1)

                ul = ui
                ur = uo
                ul(2:3) = matmul(R, ul(2:3))
                ur(2:3) = matmul(R, ur(2:3))



                fl = fEuler(ul)
                fr = fEuler(ur)

                vl = ul(2:3)/ul(1)
                vr = ur(2:3)/ur(1)

                al = min(vl(1), vr(1)) - max(anum(ul), anum(ur)) - 0.1
                ar = max(vl(1), vr(1)) + max(anum(ul), anum(ur)) + 0.1

                ! Calculate HLL intermediate state and slope
                um = (-al*ul + ar*ur + fl - fr)/(ar-al)

                ! Calculate the intercell flux
                if (0.0D0 < al) then                                                            ! Supersonic from left
                        EulerHLL = fl
                elseif (0.0D0 < ar) then                                                        ! Transonic case
                        EulerHLL = (ar*fl - al*fr + al*ar*(ur - ul)) / (ar - al)
                else                                                                            ! Supersonic from the right
                        EulerHLL = fr
                end if
                !EulerHLL(2:3) =matmul(Ri, EulerHLL(2:3))

        end function EulerHLLrot


        ! Function "mirrors" the conserved variables, i.e. the 
        ! reflection at a hard wall is mimicked by reflecting the moment while
        ! keeping density and pressure (and therefore internal energy) as is
        function MirrorConsVar(u, n) result (ur)
                double precision, dimension(4), intent(in) :: u 
                double precision, dimension(2), intent(in) :: n
                double precision, dimension(4) :: ur
                ur(:) = u(:)
                ur(2:3) = ur(2:3) - 2*n*dot_product(n, u(2:3))/dot_product(n, n)
        end function MirrorConsVar

        ! The fluxes used at a reflecting wall - x direction
        function fEulerRefl(u, n)
                double precision, dimension(4), intent(in) :: u
                double precision, dimension(2), intent(in) :: n
                double precision, dimension(4) :: fEulerRefl
                fEulerRefl = fEulerLLF(u, MirrorConsVar(u, n), n)
        end function fEulerRefl

        ! y direction
        function gEulerRefl(u, n)
                double precision, dimension(4), intent(in) :: u
                double precision, dimension(2), intent(in) :: n
                double precision, dimension(4) :: gEulerRefl
                gEulerRefl = gEulerLLF(u, MirrorConsVar(u, n), n)
        end function gEulerRefl

        ! only normal flux, LLF style
        function EulerReflLLF(u, n)
                double precision, dimension(4), intent(in) :: u
                double precision, dimension(2), intent(in) :: n
                double precision, dimension(4) :: EulerReflLLF
                EulerReflLLF = EulerLLF(u, MirrorConsVar(u, n), n)
        end function EulerReflLLF

        ! only normal flux, HLL style
        function EulerReflHLL(u, n)
                double precision, dimension(4), intent(in) :: u
                double precision, dimension(2), intent(in) :: n
                double precision, dimension(4) :: EulerReflHLL
                EulerReflHLL = EulerHLLnf(u, MirrorConsVar(u, n), n)
        end function EulerReflHLL


        ! Calculates Physical Entropy for the Euler System
        function PUEuler(u)
                implicit none
                double precision, dimension(4), intent(in) :: u
                double precision :: PUEuler
                double precision :: S
                S = log(p(u)*u(1)**(-gamma))
                PUEuler = -u(1)*S
        end function PUEuler

        ! Entropy variables
        function dPUEuler(u)
                implicit none
                double precision, dimension(4), intent(in) :: u
                double precision, dimension(4) :: dPUEuler
                double precision :: S
                double precision, dimension(4) :: dS
                S = log(p(u)*u(1)**(-gamma))
                dS = 1.0D0/(p(u)*u(1)**(-gamma))*dp(u)*u(1)**(-gamma)
                dS(1) = dS(1) + 1.0D0/(p(u)*u(1)**(-gamma))*p(u)*u(1)**(-gamma-1)*(-gamma)
                dPUEuler = -u(1)*dS
                dPUeuler(1) = dPUeuler(1) - S
        end function dPUEUler

        ! Calculates entropy flux in x direction for the physical entropy
        function PFEuler(u)
                implicit none
                double precision, dimension(4), intent(in) :: u
                double precision :: PFEuler
                double precision :: S
                S = log(p(u)*u(1)**(-gamma))
                PFEuler = -u(2)*S
        end function PFEuler

        ! Calculates Entropy flux in y direction for the physical entropy
        function PGEuler(u)
                implicit none
                double precision, dimension(4), intent(in) :: u
                double precision :: PGEuler
                double precision :: S
                S = log(p(u)*u(1)**(-gamma))
                PGEuler = -u(3)*S
        end function PGEuler

        ! Lax-Friedrichs numerical entropy flux in x direction
        function PFEulerLF(ul, ur, lambda)
                implicit none
                double precision, dimension(4), intent(in) :: ul, ur
                double precision, intent(in) :: lambda
                double precision :: PFEulerLF
                PFEulerLF = 0.5D0*(PFEuler(ul) + PFEuler(ur) - (PUEuler(ur)-PUEuler(ul))/(2*lambda))
        end function PFEulerLF

        ! Lax-Friedrichs numerical entropy flux in y direction
        function PGEulerLF(ub, ut, lambda)
                implicit none
                double precision, dimension(4), intent(in) :: ub, ut
                double precision, intent(in) :: lambda
                double precision :: PGEulerLF
                PGEulerLF = 0.5D0*(PGEuler(ub) + PGEuler(ut) - (PUEuler(ut) - PUEuler(ub))/(2*lambda))
        end function PGEulerLF

        ! Local Lax-Friedrichs numerical enropy fluxes in x
         function PFEulerLLF(ui, uo, n)
                double precision :: PFEulerLLF
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                PFEulerLLF = (PFEuler(ui) + PFEuler(uo) + 1.0*n(1)/norm2(n)*cmax(ui, uo)*(PUEuler(ui) - PUEuler(uo)))/2.0D0
        end function PFEulerLLF

        ! local Lax-Friedrichs flux in y
        function PGEulerLLF(ui, uo, n)
                double precision :: PGEulerLLF
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                PGEulerLLF = (PGEuler(ui) + PGEuler(uo) + 1.0*n(2)/norm2(n)*cmax(ui, uo)*(PUEuler(ui) - PUEuler(uo)))/2.0D0
        end function PGEulerLLF

        ! Face sided LLF entropy flux
        function PEulerLLF(ui, uo, n)
                double precision :: PEulerLLF
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                double precision, dimension(2) :: nn
                nn = n / norm2(n)
                PEulerLLF = 0.5D0*(nn(1)*(PFEuler(ui) + PFEuler(uo)) + nn(2)*(PGEuler(ui) + PGEuler(uo)) &
                        + cmax(ui, uo)*(PUEuler(ui) - PUEuler(uo)))
        end function PEulerLLF

        function PEulerHLL(ui, uo, n)
                implicit none
                double precision :: PEulerHLL
                double precision, dimension(4) :: ui, uo
                double precision, dimension(2) :: n
                ! Aux
                double precision, dimension(2) :: nn, no
                double precision, dimension(4) :: fi, fo, um
                double precision :: ai, ao, vi, vo, pfi, pfo, fis, fos
                double precision :: mso, msi
                nn = n / norm2(n)
                no(1) = -nn(2)
                no(2) = nn(1)

                fi = nn(1)*fEuler(ui) + nn(2)*gEuler(ui)
                pfi = nn(1)*PFEuler(ui) + nn(2)*PGEuler(ui)
                fo = nn(1)*fEuler(uo) + nn(2)*gEuler(uo)
                pfo = nn(1)*PFEuler(uo) + nn(2)*PGEuler(uo)
                vi = dot_product(nn, ui(2:3))/ui(1)
                vo = dot_product(nn, uo(2:3))/uo(1)

                msi = dot_product(no, ui(2:3))/ui(1)
                mso = dot_product(no, uo(2:3))/uo(1)

                ai = min(vi, vo) - max(anum(ui), anum(uo)) - abs(msi - mso)
                ao = max(vi, vo) + max(anum(ui), anum(uo)) + abs(msi - mso)


                ! Calculate HLL intermediate state and slope
                um = (-ai*ui + ao*uo + fi - fo)/(ao-ai)

                ! Calculate the intercell flux
                if (0.0D0 < ai) then                                                            ! Supersonic from left
                        PEulerHLL = pfi
                elseif (0.0D0 < ao) then                                                        ! Transonic case
                        fis = pfi + ai*(PUEuler(um) - PUEuler(ui))
                        fos = pfo + ao*(PUEuler(um) - PUEuler(uo))
                        PEulerHLL = 0.5*(fis + fos)
                else                                                                            ! Supersonic from the right
                        PEulerHLL = pfo
                end if
                
        end function PEulerHLL


        ! An entropy inequality prediction using the one-dimensional theory
        function EulerDispOpLLF(ui, uo, n) result (sigma)
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                double precision :: sigma
                
                double precision :: cb
                double precision, dimension(2) :: nn
                double precision, dimension (4) :: fi, fo, um 
                double precision :: Efi, Efo

                nn = n / norm2(n)
                cb = cmax(ui, uo)
                fi = nn(1)*fEuler(ui) + nn(2)*gEuler(ui)
                fo = nn(1)*fEuler(uo) + nn(2)*gEuler(uo)
                um = 0.5*(ui + uo) + (fi - fo) / (2*cb)
                Efi = nn(1)*PFEuler(ui) + nn(2)*PGEuler(ui)
                Efo = nn(1)*PFEuler(uo) + nn(2)*PGEuler(uo)
                sigma = cb*(2*PUEuler(um) - PUEuler(ui) - PUEuler(uo)) - (Efi - Efo)
        end function EulerDispOpLLF

        ! An entropy inequality prediction using the one-dimensional theory
        function EulerDispOpHLL(ui, uo, n) result (sigma)
                double precision, dimension(4), intent(in) :: ui, uo
                double precision, dimension(2), intent(in) :: n
                double precision :: sigma
                    ! Auxillaries
                double precision, dimension(2) :: nn, no
                double precision, dimension(4) :: fi, fo, um
                double precision :: ai, ao, vi, vo, pfi, pfo
                double precision :: msi, mso
                nn = n / norm2(n)
                no(1) = -nn(2)
                no(2) = nn(1)


                fi = nn(1)*fEuler(ui) + nn(2)*gEuler(ui)
                pfi = nn(1)*PFEuler(ui) + nn(2)*PGEuler(ui)
                fo = nn(1)*fEuler(uo) + nn(2)*gEuler(uo)
                pfo = nn(1)*PFEuler(uo) + nn(2)*PGEuler(uo)
                vi = dot_product(nn, ui(2:3))/ui(1)
                vo = dot_product(nn, uo(2:3))/uo(1)
                msi = dot_product(no, ui(2:3))/ui(1)
                mso = dot_product(no, uo(2:3))/uo(1)

                ai = min(vi, vo) - max(anum(ui), anum(uo)) - abs(msi - mso)
                ao = max(vi, vo) + max(anum(ui), anum(uo)) + abs(msi - mso)


                ! Calculate HLL intermediate state and slope
                um = (-ai*ui + ao*uo + fi - fo)/(ao-ai)

                sigma = (ao- ai)*PUEuler(um) + ai*PUEuler(ui) - ao*PUEuler(uo) - (pfi - pfo)
        end function EulerDispOpHLL

end module euler

