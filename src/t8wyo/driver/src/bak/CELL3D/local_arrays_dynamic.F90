module local_arrays_dynamic
  use my_kinddefs
  implicit none

      integer(i4)                :: mgroup
      parameter(mgroup = 100)

      integer(i4)    ntetra,npyr,nprizm,nhex,nnode,nbnode
      integer(i4)    nbface3,nbface4
      integer(i4)    nface3,nface4
      integer(i4)    ntetra_ng,npyr_ng,nprizm_ng,nhex_ng,ngpt
      integer(i4)    npatch,ncomp,nbod
      integer(i4)    nbodytag

      integer(i4) , pointer     :: ndc4(:)       => null()
      integer(i4) , pointer     :: ndc5(:)       => null()
      integer(i4) , pointer     :: ndc6(:)       => null()
      integer(i4) , pointer     :: ndc8(:)       => null()
      integer(i4) , pointer     :: nbf3(:)       => null()
      integer(i4) , pointer     :: nbf4(:)       => null()
      integer(i4) , pointer     :: ifpat3(:)     => null()
      integer(i4) , pointer     :: ifpat4(:)     => null()
      real(r8)    , pointer     :: xgeom(:)      => null()
      integer(i4) , pointer     :: ndf3(:)       => null()
      integer(i4) , pointer     :: ndf4(:)       => null()

      CHARACTER*20, pointer  :: ccomp(:)          => null()
      CHARACTER*20, pointer  :: cbod(:)           => null()

      integer(i4)  , pointer     :: nde(:)        => null()
      real(r8)     , pointer     :: gnorm(:)      => null()
      real(r8)     , pointer     :: vcoef(:)      => null()
      real(r8)     , pointer     :: vol(:)        => null()
      real(r8)     , pointer     :: eps(:)        => null()
      real(r8)     , pointer     :: feps(:)       => null()

      integer(i4)  , pointer     :: ibnode(:)     => null()
      integer(i4)  , pointer     :: ibinv(:)      => null()
      real(r8)     , pointer     :: fbcnorm(:)    => null()
      integer(i4)  , pointer     :: kpatch_2_bcs(:) => null()

      integer(i4)  , pointer     :: ilvertex(:)   => null()
      integer(i4)  , pointer     :: jll(:)        => null()
      integer(i4)  , pointer     :: nlvec(:)      => null()
      integer(i4)                :: klgroup

      integer(i4)                :: ndegrp 
      integer(i4)  , pointer     :: ndevec(:)     => null()

      integer(i4)  , pointer     :: iinodev(:)    => null()
      integer(i4)  , pointer     :: jjnodev(:)    => null()
      integer(i4)  , pointer     :: iwgt(:)       => null()

      integer(i4)  , pointer     :: kgpt_id(:)    => null()
      integer(i4)  , pointer     :: kgpt_ip(:)    => null()


end module local_arrays_dynamic
