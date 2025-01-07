! File contains the fluxes for the two-dimensional euler equations and the corresponding entropy fluxes
module euler
        use param
        implicit none
        real(np) :: gamma = 1.4_np
        logical, parameter :: kse = .true.     ! Switches to more complex speeds
        contains

        ! calculates pressure from conserved variables
        function p(u)
                implicit none
                real(np) :: u(4)
                real(np) :: p
                p = (gamma-1.0_np)*(u(4)-0.5_np*(u(2)**2 + u(3)**2)/u(1))
        end function p

        ! derivative of the pressure from conserved variables
        function dp(u)
                implicit none
                real(np) :: u(4)
                real(np) :: dp(4)
                dp(1) = (gamma-1)*0.5_np*(u(2)**2 + u(3)**2)*u(1)**(-2.0_np)
                dp(2) = -(gamma-1)*(u(2)/u(1))
                dp(3) = -(gamma-1)*(u(3)/u(1))
                dp(4) = (gamma-1)
        end function dp

        ! calculates speed of sound from conserved variables
        function a(u)
                implicit none
                real(np) :: u(4)
                real(np) :: a
                a = sqrt(gamma*p(u)/ u(1))
        end function a

        ! Calcualtes the 'numeric' sound speed - used to remove problems with
        ! sonic point glitches following Tang:"On the sonic point glitch"
        ! in the HLL numerical flux.
        function anum(u)
                implicit none
                real(np) :: u(4)
                real(np) :: anum
                anum = sqrt(gamma*p(u)/ u(1))
        end function anum


        ! speed of sound directly from pressure and density
        function apr(press, rho)
                implicit none
                real(np) :: press, rho
                real(np) :: apr
                apr = sqrt(gamma*press/rho)
        end function apr

        ! 'numerical' speed of sound directly from pressure and density
        function aprnum(press, rho)
                implicit none
                real(np) :: press, rho
                real(np) :: aprnum
                aprnum = sqrt(gamma*press/rho)
        end function aprnum

        ! Converts between physical and conserved variables
        function ConsVar(rho, vx, vy, p)
                implicit none
                real(np) :: ConsVar(4)
                real(np) :: rho, vx, vy, p
                ConsVar(1) = rho
                ConsVar(2) = rho*vx
                ConsVar(3) = rho*vy
                ConsVar(4) = p / (gamma-1.0_np) + 0.5_np * (vx**2 + vy**2)*rho
        end function ConsVar

        ! Converts from conserved variables to physical variables
        function PhysVar(rho, rhovx, rhovy, E)
                implicit none
                real(np) :: PhysVar(4)
                real(np) :: rho, rhovx, rhovy, E
                PhysVar(1) = rho
                PhysVar(2) = rhovx/rho
                PhysVar(3) = rhovy/rho
                PhysVar(4) = (gamma-1.0_np)*(E - 0.5_np * (rhovx**2 + rhovy**2) / rho)
        end function PhysVar


        ! The following function relates to Section 3.1.3 of
        ! Toro: Riemann Solvers and Numericl Methods.
        ! given right data it calculates the left states for a
        ! shock moving with Mach number Ms
        !       rhor: right density
        !       vr: right speed
        !       pr: right pressure
        !       Ms: Mach number of the shock
        ! returns u in conserved variables 
        function Lvars(rhor, vr, pr, Ms) result (u)
                real(np) :: rhor, vr, pr, Ms
                real(np) :: u(4)
                ! Aux
                real(np) :: Mr, rhol, pl, vl, S
                ! Calculate the right mach number
                Mr = vr / sqrt(gamma*pr/rhor)
                ! Calculate the Shock speed concerning the right sound speed
                S = Ms * sqrt(gamma*pr/rhor)
                rhol = rhor*(gamma + 1.0_np)*(Mr - Ms)**2 &
                        / ((gamma - 1.0_np)*(Mr - Ms)**2 + 2.0_np)
                pl = pr * (2.0_np * gamma*(Mr - Ms)**2 - gamma + 1.0_np) &
                        / (gamma + 1.0_np)
                vl = (1.0_np - rhor/rhol)*S + vr*rhor/rhol
                u = ConsVar(rhol, vl, 0.0_np, pl)
        end function Lvars

        ! calculates x component of the Euler flux
        function fEuler(u)
                implicit none
                real(np) :: u(4)
                real(np) :: fEuler(4)
                real(np) :: pu
                pu = p(u)
                fEuler(1) = u(2)
                fEuler(2) = u(2)**2/u(1) + pu
                fEuler(3) = u(2)*u(3)/u(1)
                fEuler(4) = u(2)/u(1)*(u(4) + pu)
        end function fEuler

        ! calculates y component of Euler flux
        function gEuler(u)
                implicit none
                real(np) :: u(4)
                real(np) :: gEuler(4)
                real(np) :: pu
                pu = p(u)
                gEuler(1) = u(3)
                gEuler(2) = u(2)*u(3)/u(1)
                gEuler(3) = u(3)**2/u(1) + pu
                gEuler(4) = u(3)/u(1)*(u(4) + pu)
        end function gEuler

        ! calculates Lax-Friedrichs flux in x
        function fEulerLF(ui, uo, n, lambda)
                real(np) :: fEulerLF(4)
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                real(np), intent(in) :: lambda
                fEulerLF = (fEuler(ui) + fEuler(uo) + sign(1.0_np, n(1))*(ui - uo)/(2.0_np*lambda))/2.0_np
        end function fEulerLF

        ! Lax-Friedrichs flux in y
        function gEulerLF(ui, uo, n, lambda)
                real(np) :: gEulerLF(4)
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                real(np), intent(in) :: lambda
                gEulerLF = (gEuler(ui) + gEuler(uo) + sign(1.0_np, n(2))*(ui - uo)/(2*lambda))/2.0_np
        end function gEulerLF

        ! calculates the maximum speed for the conection of two given states
        function cmax(ul, ur)
                real(np) :: ul(4), ur(4)
                real(np) :: cmax
                real(np) :: vl(2), vr(2)
                vl = ul(2:3) / ul(1)
                vr = ur(2:3) / ur(1)
                cmax = max(norm2(vl), norm2(vr)) + max(anum(ul), anum(ur))
        end function cmax


        ! Calcualtes local Lax-Friedrichs flux in x
        function fEulerLLF(ui, uo, n)
                real(np) :: fEulerLLF(4)
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                fEulerLLF = (fEuler(ui) + fEuler(uo) + 1.0_np*n(1)/norm2(n)*cmax(ui, uo)*(ui - uo))/2.0_np
        end function fEulerLLF

        ! local Lax-Friedrichs flux in y
        function gEulerLLF(ui, uo, n)
                real(np) :: gEulerLLF(4)
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                gEulerLLF = (gEuler(ui) + gEuler(uo) + 1.0_np*n(2)/norm2(n)*cmax(ui, uo)*(ui - uo))/2.0_np
        end function gEulerLLF

        ! Normal Lax-Friedrichs flux
        function EulerLLF(ui, uo, n)
                real(np) :: EulerLLF(4)
                real(np) :: ui(4), uo(4)
                real(np) :: n(2)
                real(np) :: nn(2)
                nn = n / norm2(n)
                EulerLLF = 0.5_np*(nn(1)*(fEuler(ui) + fEuler(uo)) + nn(2)*(gEuler(ui) + gEuler(uo)) &
                        + cmax(ui, uo)*(ui - uo))
        end function EulerLLF

        ! direct computation of the HLL flux
        function EulerHLLnf(ui, uo, n) result (EulerHLL)
                implicit none
                real(np) :: EulerHLL(4)
                real(np) :: ui(4), uo(4)
                real(np) :: n(2)
                ! Auxillaries
                real(np) :: nn(2), no(2)
                real(np) :: fi(4), fo(4)
                real(np) :: ai, ao, vi, vo, msi, mso
                real(np) :: anumb, pi, po, vii, vio, voi, voo
                nn = n / norm2(n)
                no(1) = -nn(2)
                no(2) = nn(1)

                anumb = aprnum(max(p(ui), p(uo)), min(ui(1), uo(1)))

                fi = nn(1)*fEuler(ui) + nn(2)*gEuler(ui)
                fo = nn(1)*fEuler(uo) + nn(2)*gEuler(uo)
                msi = dot_product(no, ui(2:3)/ui(1))
                mso = dot_product(no, uo(2:3)/uo(1))

                if (kse) then
                        pi = dot_product(nn, ui(2:3))
                        po = dot_product(nn, uo(2:3))
                        vii = pi / ui(1)
                        vio = pi / uo(1)
                        voi = po / ui(1)
                        voo = po / uo(1)

                        ai = min(vii, vio, voi, voo) - anumb - 0.5_np*abs(msi - mso)
                        ao = max(vii, vio, voi, voo) + anumb + 0.5_np*abs(msi - mso)

                else
                        vi = dot_product(nn, ui(2:3))/ui(1)
                        vo = dot_product(nn, uo(2:3))/uo(1)

                        ai = min(vi, vo) - max(anum(ui), anum(uo)) - 0.5_np*abs(msi - mso)
                        ao = max(vi, vo) + max(anum(ui), anum(uo)) + 0.5_np*abs(msi - mso)
                end if


                if (0.0_np < ai) then                                                            ! Supersonic from left
                        EulerHLL = fi
                elseif (0.0_np < ao) then                                                        ! Transonic case
                        EulerHLL = (ao*fi - ai*fo + ai*ao*(uo - ui)) / (ao - ai)
                else                                                                            ! Supersonic from the right
                        EulerHLL = fo
                end if
        end function EulerHLLnf

        ! Computation using the rotation onto the x direction
        function EulerHLLrot(ui, uo, n) result (EulerHLL)
                implicit none
                real(np), dimension(4) :: EulerHLL(4)
                real(np), dimension(4) :: ui(4), uo(4)
                real(np), dimension(2) :: n(2)
                ! Auxillaries
                real(np), dimension(2) :: nn(2), vl(2), vr(2)
                real(np), dimension(4) :: fl(4), fr(4), ul(4), ur, um
                real(np), dimension(2, 2) :: R, Ri
                real(np) :: al, ar
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

                al = min(vl(1), vr(1)) - max(anum(ul), anum(ur)) - 0.1_np
                ar = max(vl(1), vr(1)) + max(anum(ul), anum(ur)) + 0.1_np

                ! Calculate HLL intermediate state and slope
                um = (-al*ul + ar*ur + fl - fr)/(ar-al)

                ! Calculate the intercell flux
                if (0.0_np < al) then                                                            ! Supersonic from left
                        EulerHLL = fl
                elseif (0.0_np < ar) then                                                        ! Transonic case
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
                real(np), intent(in) :: u(4)
                real(np), intent(in) :: n(2)
                real(np) :: ur(4)
                ur(:) = u(:)
                ur(2:3) = ur(2:3) - 2.0_np*n*dot_product(n, u(2:3))/dot_product(n, n)
        end function MirrorConsVar

        ! The fluxes used at a reflecting wall - x direction
        function fEulerRefl(u, n)
                real(np), intent(in) :: u(4)
                real(np), intent(in) :: n(2)
                real(np) :: fEulerRefl(4)
                fEulerRefl = fEulerLLF(u, MirrorConsVar(u, n), n)
        end function fEulerRefl

        ! y direction
        function gEulerRefl(u, n)
                real(np), intent(in) :: u(4)
                real(np), intent(in) :: n(2)
                real(np) :: gEulerRefl(4)
                gEulerRefl = gEulerLLF(u, MirrorConsVar(u, n), n)
        end function gEulerRefl

        ! only normal flux, LLF style
        function EulerReflLLF(u, n)
                real(np), intent(in) :: u(4)
                real(np), intent(in) :: n(2)
                real(np) :: EulerReflLLF(4)
                EulerReflLLF = EulerLLF(u, MirrorConsVar(u, n), n)
        end function EulerReflLLF

        ! only normal flux, HLL style
        function EulerReflHLL(u, n)
                real(np), intent(in) :: u(4)
                real(np), intent(in) :: n(2)
                real(np) :: EulerReflHLL(4)
                EulerReflHLL = EulerHLLnf(u, MirrorConsVar(u, n), n)
        end function EulerReflHLL

        function h(S)
                real(np) :: S, h
                h = (gamma + 1.0_np) / (gamma - 1.0_np) *exp(S / (gamma + 1.0_np))
        end function h

        function dh(S)
                real(np) :: S, dh
                dh = exp(S / (gamma + 1.0_np)) / (gamma - 1.0_np)
        end function dh


        ! Calculates the Physical Entropy for the Euler System
        function PUEuler(u)
                implicit none
                real(np), intent(in) :: u(4)
                real(np) :: PUEuler
                real(np) :: S
                S = log(p(u)*u(1)**(-gamma))
                PUEuler = -u(1)*S
        end function PUEuler

        ! Entropy variables for the physical entropy
        function dPUEuler(u)
                implicit none
                real(np), intent(in) :: u(4)
                real(np) :: dPUEuler(4)
                real(np) :: S
                real(np) :: dS(4)
                S = log(p(u)*u(1)**(-gamma))
                dS = 1.0_np/(p(u)*u(1)**(-gamma))*dp(u)*u(1)**(-gamma)
                dS(1) = dS(1) + 1.0_np/(p(u)*u(1)**(-gamma))*p(u)*u(1)**(-gamma-1)*(-gamma)
                dPUEuler = -u(1)*dS
                dPUeuler(1) = dPUeuler(1) - S
        end function dPUEUler

        ! Calculates physical entropy flux in x direction for the physical entropy
        function PFEuler(u)
                implicit none
                real(np), intent(in) :: u(4)
                real(np) :: PFEuler
                real(np) :: S
                S = log(p(u)*u(1)**(-gamma))
                PFEuler = -u(2)*S
        end function PFEuler

        ! Calculates Physical Entropy flux in y direction for the physical entropy
        function PGEuler(u)
                implicit none
                real(np), intent(in) :: u(4)
                real(np) :: PGEuler
                real(np) :: S
                S = log(p(u)*u(1)**(-gamma))
                PGEuler = -u(3)*S
        end function PGEuler

        ! Lax-Friedrichs numerical entropy flux in x direction
        function PFEulerLF(ul, ur, lambda)
                implicit none
                real(np), intent(in) :: ul(4), ur(4)
                real(np), intent(in) :: lambda
                real(np) :: PFEulerLF
                PFEulerLF = 0.5_np*(PFEuler(ul) + PFEuler(ur) - (PUEuler(ur)-PUEuler(ul))/(2.0_np*lambda))
        end function PFEulerLF

        ! Lax-Friedrichs numerical entropy flux in y direction
        function PGEulerLF(ub, ut, lambda)
                implicit none
                real(np), intent(in) :: ub(4), ut(4)
                real(np), intent(in) :: lambda
                real(np) :: PGEulerLF
                PGEulerLF = 0.5_np*(PGEuler(ub) + PGEuler(ut) - (PUEuler(ut) - PUEuler(ub))/(2.0_np*lambda))
        end function PGEulerLF

        ! Local Lax-Friedrichs numerical enropy fluxes in x
         function PFEulerLLF(ui, uo, n)
                real(np) :: PFEulerLLF
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                PFEulerLLF = (PFEuler(ui) + PFEuler(uo) + n(1)/norm2(n)*cmax(ui, uo)*(PUEuler(ui) - PUEuler(uo)))/2.0_np
        end function PFEulerLLF

        ! local Lax-Friedrichs flux in y
        function PGEulerLLF(ui, uo, n)
                real(np) :: PGEulerLLF
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                PGEulerLLF = (PGEuler(ui) + PGEuler(uo) + n(2)/norm2(n)*cmax(ui, uo)*(PUEuler(ui) - PUEuler(uo)))/2.0_np
        end function PGEulerLLF

        ! Face sided LLF entropy flux
        function PEulerLLF(ui, uo, n)
                real(np) :: PEulerLLF
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                real(np) :: nn(2)
                nn = n / norm2(n)
                PEulerLLF = 0.5_np*(nn(1)*(PFEuler(ui) + PFEuler(uo)) + nn(2)*(PGEuler(ui) + PGEuler(uo)) &
                        + cmax(ui, uo)*(PUEuler(ui) - PUEuler(uo)))
        end function PEulerLLF

        ! HLL entropy flux for the physical entropy
        function PEulerHLL(ui, uo, n)
                implicit none
                real(np) :: PEulerHLL
                real(np) :: ui(4), uo(4)
                real(np) :: n(2)
                ! Aux
                real(np) :: nn(2), no(2)
                real(np) :: fi(4), fo(4), um(4)
                real(np) :: ai, ao, vi, vo, pfi, pfo, fis, fos
                real(np) :: mso, msi
                real(np) :: anumb, pi, po, vii, vio, voi, voo
                anumb = aprnum(max(p(ui), p(uo)), min(ui(1), uo(1)))
                nn = n / norm2(n)
                no(1) = -nn(2)
                no(2) = nn(1)

                fi = nn(1)*fEuler(ui) + nn(2)*gEuler(ui)
                pfi = nn(1)*PFEuler(ui) + nn(2)*PGEuler(ui)
                fo = nn(1)*fEuler(uo) + nn(2)*gEuler(uo)
                pfo = nn(1)*PFEuler(uo) + nn(2)*PGEuler(uo)

                
                msi = dot_product(no, ui(2:3))/ui(1)
                mso = dot_product(no, uo(2:3))/uo(1)

                if (kse) then
                        pi = dot_product(nn, ui(2:3))
                        po = dot_product(nn, uo(2:3))
                        vii = pi / ui(1)
                        vio = pi / uo(1)
                        voi = po / ui(1)
                        voo = po / uo(1)

                        ai = min(vii, vio, voi, voo) - anumb - 0.5_np*abs(msi - mso)
                        ao = max(vii, vio, voi, voo) + anumb + 0.5_np*abs(msi - mso)
                else
                        vi = dot_product(nn, ui(2:3))/ui(1)
                        vo = dot_product(nn, uo(2:3))/uo(1)

                        ai = min(vi, vo) - max(anum(ui), anum(uo)) - 0.5_np*abs(msi - mso)
                        ao = max(vi, vo) + max(anum(ui), anum(uo)) + 0.5_np*abs(msi - mso)
                end if

                ! Calculate HLL intermediate state and slope
                um = (-ai*ui + ao*uo + fi - fo)/(ao-ai)

                ! Calculate the intercell flux
                if (0.0_np < ai) then                                                            ! Supersonic from left
                        PEulerHLL = pfi
                elseif (0.0_np < ao) then                                                        ! Transonic case
                        fis = pfi + ai*(PUEuler(um) - PUEuler(ui))
                        fos = pfo + ao*(PUEuler(um) - PUEuler(uo))
                        PEulerHLL = 0.5_np*(fis + fos)
                else                                                                            ! Supersonic from the right
                        PEulerHLL = pfo
                end if
                
        end function PEulerHLL

        ! An entropy inequality prediction using the one-dimensional theory
        function EulerDispOpLLF(ui, uo, n) result (sigma)
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                real(np) :: sigma
                
                real(np) :: cb
                real(np) :: nn(2)
                real(np) :: fi(4), fo(4), um(4) 
                real(np) :: Efi, Efo

                nn = n / norm2(n)
                cb = cmax(ui, uo)
                fi = nn(1)*fEuler(ui) + nn(2)*gEuler(ui)
                fo = nn(1)*fEuler(uo) + nn(2)*gEuler(uo)
                um = 0.5_np*(ui + uo) + (fi - fo) / (2*cb)
                Efi = nn(1)*PFEuler(ui) + nn(2)*PGEuler(ui)
                Efo = nn(1)*PFEuler(uo) + nn(2)*PGEuler(uo)
                sigma = cb*(2*PUEuler(um) - PUEuler(ui) - PUEuler(uo)) - (Efi - Efo)
        end function EulerDispOpLLF

        ! An entropy inequality prediction using the one-dimensional theory
        function EulerDispOpHLL(ui, uo, n) result (sigma)
                real(np), intent(in) :: ui(4), uo(4)
                real(np), intent(in) :: n(2)
                real(np) :: sigma
                    ! Auxillaries
                real(np) :: nn(2), no(2)
                real(np) :: fi(4), fo(4), um(4)
                real(np) :: ai, ao, vi, vo, pfi, pfo
                real(np) :: msi, mso
                real(np) :: anumb, pi, po, vii, vio, voi, voo
                nn = n / norm2(n)
                no(1) = -nn(2)
                no(2) = nn(1)

                anumb = aprnum(max(p(ui), p(uo)), min(ui(1), uo(1)))

                fi = nn(1)*fEuler(ui) + nn(2)*gEuler(ui)
                pfi = nn(1)*PFEuler(ui) + nn(2)*PGEuler(ui)
                fo = nn(1)*fEuler(uo) + nn(2)*gEuler(uo)
                pfo = nn(1)*PFEuler(uo) + nn(2)*PGEuler(uo)
                msi = dot_product(no, ui(2:3))/ui(1)
                mso = dot_product(no, uo(2:3))/uo(1)

                if (kse) then
                        pi = dot_product(nn, ui(2:3))
                        po = dot_product(nn, uo(2:3))
                        vii = pi / ui(1)
                        vio = pi / uo(1)
                        voi = po / ui(1)
                        voo = po / uo(1)

                        ai = min(vii, vio, voi, voo) - anumb - 0.5_np*abs(msi - mso)
                        ao = max(vii, vio, voi, voo) + anumb + 0.5_np*abs(msi - mso)
                else
                        vi = dot_product(nn, ui(2:3))/ui(1)
                        vo = dot_product(nn, uo(2:3))/uo(1)

                        ai = min(vi, vo) - max(anum(ui), anum(uo)) - 0.5_np*abs(msi - mso)
                        ao = max(vi, vo) + max(anum(ui), anum(uo)) + 0.5_np*abs(msi - mso)
                end if

                ! Calculate HLL intermediate state
                um = (-ai*ui + ao*uo + fi - fo)/(ao-ai)

                sigma = (ao- ai)*PUEuler(um) + ai*PUEuler(ui) - ao*PUEuler(uo) - (pfi - pfo)
        end function EulerDispOpHLL

end module euler

