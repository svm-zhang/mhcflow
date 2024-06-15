`polysolvermod` uses [mamba](https://github.com/mamba-org/mamba) and [boa](https://github.com/mamba-org/boa) as package manager and builder.

## Linux-x64

To install `polysolvermod` on a LINUX platform, simply run the following:

```
cd "$polysolvermod_repo"
boa build .     # this builds polysolvermod as a local tarball
mamba create -n hlatyping
mamba install -n hlatyping --use-local polysolvermod
```

One important thing to note: you will need to have `novoalign` and `novoindex` binaries available on your `PATH`. `polysolvermod` does not provide them and you can download them from novocraft website.

## OSX-64 (Apple Silicon chip)

Installation on OSX-64 requires a few additional steps:
* You need to switch to `bash` shell
* Set `subdirs` to `osx-64` in your `.condarc` or `.mambarc` config file
* You need install GNU-coreutils to have access to the GNU implementation of `grep`, `split`, `cut`, `date`, and `du`
```
brew install coreutils
```

`homebrew` adds a `g` prefix to these commands to avoid conflict with the Free-BSD ones on Mac. To make sure `polysolvermod` to use the GNU version, you can simply

```
cp "$homebrew/bin/ggrep" "$conda_env_dir/bin/grep"
cp "$homebrew/bin/gcut" "$conda_env_dir/bin/cut"
cp "$homebrew/bin/gsplit" "$conda_env_dir/bin/split"
cp "$homebrew/bin/gdate" "$conda_env_dir/bin/date"
cp "$homebrew/bin/gdu" "$conda_env_dir/bin/du"
```

You can also use symbolic link rather than copying.

One more thing: according to novocraft website, `novoalign` are yet to be optimized for the Apple Silicon chip.

## Manual

I try to minimize number of dependencies required to run `polysolvermod`. Part of the reason I decide to completely re-write typer in python is to remove the `perl` dependency that causes me so much trouble when installing the original `polysolver` (let me know if you want to have a `recipe.yaml` for `polysolver`).

Please take a look at the runtime depencies defined in the `recipe.yaml` file in the repo, and have them installed on your local environment. Next, you install `polysolvermod` by simply running a few commands:

```
cd "$polysolvermod_repo"
mkdir bin
cp scripts/hlapolysolver.sh bin/polysolvermod
cp scripts/polysolver_realigner.sh bin/realigner
cp scripts/fisher.sh bin/fisher
cp scripts/extract_sample_hlaref.sh bin/extractor
cp "$novocraft_dir/novoalign" bin/novoalign     # you dont need this if you have novoalign on your PATH
cp "$novocraft_dir/novoindex" bin/novoindex
export PATH="$polysolvermod_repo/bin:$PATH"
```

To have access to all the binaries every time you initiate a shell, you can add the last command above to your `.bash_profile` file.


## Limitations

* I'd like to keep `polysolvermod` as an offline package for now. Once people think it is helpful, I will publish it to the public conda domain
* I am aware of new successor package manager and builder to `mamba` and `boa`: [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) and [rattler](https://github.com/prefix-dev/rattler-build). As `--use-local` is [yet to be supported](https://github.com/mamba-org/mamba/issues/1991), I will stay with `mamba` for now. If anything changes, I will consider to update it in the future
* I have only tested the installation on Linux-64 and OSX-64. No WINDOWS OS supported.
