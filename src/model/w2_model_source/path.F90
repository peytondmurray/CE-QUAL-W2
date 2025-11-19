module path
    implicit none

    interface
        integer function c_chdir(path) bind(C,name="chdir")
            use iso_c_binding
            character(kind=c_char) :: path(*)
        end function
    end interface

contains
    function join_path(root, path) result(joined)
        implicit none
        character(*), intent(in)  :: root, path
        character(:), allocatable :: joined, sep
        character(:), allocatable :: cleaned_root, cleaned_path
        integer :: len_root

        cleaned_root = trim(adjustl(root))
        cleaned_path = trim(adjustl(path))

        len_root = len(cleaned_root)

#ifdef _WIN32
        sep = '\'
#else
        sep = '/'
#endif

        if (cleaned_root(len_root:len_root) == sep) then
            cleaned_root = cleaned_root(1:len_root-1)
        endif

        joined = cleaned_root//sep//cleaned_path
    end function join_path

    integer function copy_file(src, dst) result(stat)
      implicit none
      character(*), intent(in) :: src, dst
      character(len=512) :: msg
      character(:), allocatable :: cmd

#ifdef _WIN32
      cmd = 'copy '
#else
      cmd = 'cp '
#endif

      call execute_command_line(cmd//trim(adjustl(src))//'" "'//trim(adjustl(dst))//'"', exitstat=stat, cmdmsg=msg)

      if (stat /= 0) then
        write(*,*) 'Error writing file: ', trim(msg)
        write(9911,*) 'Error writing file: ', trim(msg)
      endif
    end function copy_file

    ! Fortran chdir. Invokes the chdir POSIX call.
    ! https://stackoverflow.com/a/26731789/8100451
    integer function f_chdir(path) result(err)
        use iso_c_binding
        character(*), intent(in) :: path

        err = c_chdir(path//c_null_char)
    end function f_chdir

    ! Make a directory using mkdir.
    integer function f_mkdir(path) result(stat)
        character(*), intent(in) :: path
        call execute_command_line("mkdir -p " // path, wait=.true., exitstat=stat)
    end function f_mkdir
end module path
