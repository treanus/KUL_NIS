# The installer moved

`setup_environment.sh`, `SOFTWARE_ROOT_SETUP.md` and `MIGRATION_GUIDE.md` used to
live here. They now have their own repository:

**https://github.com/Rad-dude/KUL_Linux_setup**

```bash
git clone https://github.com/Rad-dude/KUL_Linux_setup
cd KUL_Linux_setup
./setup_environment.sh
```

## Why

The installer clones KUL_NIS itself (along with KUL_VBG and KUL_FWT, at pinned
branches). Living inside the repo it installs meant you had to clone KUL_NIS by
hand first just to reach the script — a bootstrap inversion. It also means the
installer can now be versioned independently: pin bumps no longer land as commits
in a pipeline repo.

## Where the history is

The files were split out unchanged at KUL_NIS commit **`e9bcf40`**. Anything
before that point is in this repository's history under
`tools/setup_scripts_new/`; anything after is in KUL_Linux_setup.
