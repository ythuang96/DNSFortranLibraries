module matrices
    use types, only: sp, dp
    implicit none
    private
    include 'parameters'

    real(kind=dp),public :: dt11(5,my),dt12(5,my)

end module matrices
