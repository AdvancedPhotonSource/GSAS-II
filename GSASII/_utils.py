# ### This is copied from the Python 3.10 source
# Lib/distutils/file_util.py
# PSF license


import os


class DistutilsFileError(Exception): ...

# for generating verbose output in 'copy_file()'
_copy_action = { None:   'copying',
                 'hard': 'hard linking',
                 'sym':  'symbolically linking' }


def _copy_file(src, dst, preserve_mode=1, preserve_times=1, update=0,
              link=None, verbose=1, dry_run=0):
    """Copy a file 'src' to 'dst'.  If 'dst' is a directory, then 'src' is
    copied there with the same name; otherwise, it must be a filename.  (If
    the file exists, it will be ruthlessly clobbered.)  If 'preserve_mode'
    is true (the default), the file's mode (type and permission bits, or
    whatever is analogous on the current platform) is copied.  If
    'preserve_times' is true (the default), the last-modified and
    last-access times are copied as well.  If 'update' is true, 'src' will
    only be copied if 'dst' does not exist, or if 'dst' does exist but is
    older than 'src'.

    'link' allows you to make hard links (os.link) or symbolic links
    (os.symlink) instead of copying: set it to "hard" or "sym"; if it is
    None (the default), files are copied.  Don't set 'link' on systems that
    don't support it: 'copy_file()' doesn't check if hard or symbolic
    linking is available. If hardlink fails, falls back to
    _copy_file_contents().

    Under Mac OS, uses the native file copy function in macostools; on
    other systems, uses '_copy_file_contents()' to copy file contents.

    Return a tuple (dest_name, copied): 'dest_name' is the actual name of
    the output file, and 'copied' is true if the file was copied (or would
    have been copied, if 'dry_run' true).
    """
    # XXX if the destination file already exists, we clobber it if
    # copying, but blow up if linking.  Hmmm.  And I don't know what
    # macostools.copyfile() does.  Should definitely be consistent, and
    # should probably blow up if destination exists and we would be
    # changing it (ie. it's not already a hard/soft link to src OR
    # (not update) and (src newer than dst).

    # replaces an import in original source
    def newer (source, target):
        """Return true if 'source' exists and is more recently modified than
        'target', or if 'source' exists and 'target' doesn't.  Return false if
        both exist and 'target' is the same age or younger than 'source'.
        Raise DistutilsFileError if 'source' does not exist.
        """
        if not os.path.exists(source):
            raise DistutilsFileError("file '%s' does not exist" %
                                     os.path.abspath(source))
        if not os.path.exists(target):
            return 1

        from stat import ST_MTIME
        mtime1 = os.stat(source)[ST_MTIME]
        mtime2 = os.stat(target)[ST_MTIME]

        return mtime1 > mtime2

    from stat import ST_ATIME, ST_MTIME, ST_MODE, S_IMODE

    if not os.path.isfile(src):
        raise DistutilsFileError(
              "can't copy '%s': doesn't exist or not a regular file" % src)

    if os.path.isdir(dst):
        dir = dst
        dst = os.path.join(dst, os.path.basename(src))
    else:
        dir = os.path.dirname(dst)

    try:
        action = _copy_action[link]
    except KeyError:
        raise ValueError("invalid value '%s' for 'link' argument" % link)


    if dry_run:
        return (dst, 1)

    # If linking (hard or symbolic), use the appropriate system call
    # (Unix only, of course, but that's the caller's responsibility)
    elif link == 'hard':
        if not (os.path.exists(dst) and os.path.samefile(src, dst)):
            try:
                os.link(src, dst)
                return (dst, 1)
            except OSError:
                # If hard linking fails, fall back on copying file
                # (some special filesystems don't support hard linking
                #  even under Unix, see issue #8876).
                pass
    elif link == 'sym':
        if not (os.path.exists(dst) and os.path.samefile(src, dst)):
            os.symlink(src, dst)
            return (dst, 1)

    class DistutilsFileError(Exception):
        ...

    # otherwise (non-Mac, not linking), copy the file contents and
    # (optionally) copy the times and mode.
    def _copy_file_contents(src, dst, buffer_size=16*1024):
        """Copy the file 'src' to 'dst'; both must be filenames.  Any error
        opening either file, reading from 'src', or writing to 'dst', raises
        DistutilsFileError.  Data is read/written in chunks of 'buffer_size'
        bytes (default 16k).  No attempt is made to handle anything apart from
        regular files.
        """
        # Stolen from shutil module in the standard library, but with
        # custom error-handling added.
        fsrc = None
        fdst = None
        try:
            try:
                fsrc = open(src, 'rb')
            except OSError as e:
                raise DistutilsFileError("could not open '%s': %s" % (src, e.strerror))

            if os.path.exists(dst):
                try:
                    os.unlink(dst)
                except OSError as e:
                    raise DistutilsFileError(
                          "could not delete '%s': %s" % (dst, e.strerror))

            try:
                fdst = open(dst, 'wb')
            except OSError as e:
                raise DistutilsFileError(
                      "could not create '%s': %s" % (dst, e.strerror))

            while True:
                try:
                    buf = fsrc.read(buffer_size)
                except OSError as e:
                    raise DistutilsFileError(
                          "could not read from '%s': %s" % (src, e.strerror))

                if not buf:
                    break

                try:
                    fdst.write(buf)
                except OSError as e:
                    raise DistutilsFileError(
                          "could not write to '%s': %s" % (dst, e.strerror))
        finally:
            if fdst:
                fdst.close()
            if fsrc:
                fsrc.close()

    _copy_file_contents(src, dst)
    if preserve_mode or preserve_times:
        st = os.stat(src)

        # According to David Ascher <da@ski.org>, utime() should be done
        # before chmod() (at least under NT).
        if preserve_times:
            os.utime(dst, (st[ST_ATIME], st[ST_MTIME]))
        if preserve_mode:
            os.chmod(dst, S_IMODE(st[ST_MODE]))

    return (dst, 1)
