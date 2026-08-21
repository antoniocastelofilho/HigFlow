# Running HigFlow on Windows

HigFlow runs on Windows through WSL2 - a real Linux kernel, not an emulation
layer. Everything below was done on Windows 11 and the timings are measured, not
estimated.

```powershell
wsl --install -d Ubuntu-22.04
```

If WSL is already present on your machine, that one command is the whole
installation and needs no restart. The rest of this page is about the parts that
are not obvious.

---

## Contents

- [Why there is no native Windows build](#why-there-is-no-native-windows-build)
- [Which path to take](#which-path-to-take)
- [Setting up WSL2](#setting-up-wsl2)
- [Two things the first-run setup gets wrong](#two-things-the-first-run-setup-gets-wrong)
- [Where to put the repository](#where-to-put-the-repository)
- [Path A - Docker, for running](#path-a--docker-for-running)
- [Path B - native build in WSL, for developing](#path-b--native-build-in-wsl-for-developing)
- [Viewing results in ParaView on Windows](#viewing-results-in-paraview-on-windows)
- [Editing code from Windows](#editing-code-from-windows)
- [Disk and memory](#disk-and-memory)
- [Troubleshooting](#troubleshooting)
- [MSYS2 and Cygwin](#msys2-and-cygwin)

---

## Why there is no native Windows build

This is a limitation with specific causes, not an omission:

| Dependency | Windows |
|---|---|
| **OpenMPI** | No supported Windows build. MS-MPI exists but is a different implementation with a different ABI |
| **libfyaml** | No Windows port |
| **libnuma** | Linux-only, and `CMakeLists.txt` marks it `REQUIRED`, so configuration fails before a single file compiles |
| Zoltan (Trilinos) | Buildable, laborious |
| PETSc | Buildable through MSYS2, laborious |

The code itself is close to portable - three POSIX-only headers across 176,000
lines, no `fork`, no nested functions. It is the dependency stack that does not
cross over.

WSL2 is not a workaround for this. It is a Linux kernel running on your machine,
so the Linux build is the real build.

## Which path to take

| | Path A - Docker | Path B - native in WSL |
|---|---|---|
| Just want to run simulations | ✅ | |
| Want to modify the solver | | ✅ |
| Setup effort | one build, ~7 min | dependencies, ~30 min |
| Rebuild after editing a `.c` | rebuild image, ~40 s | `make`, seconds |
| Debugger | needs the development image | works directly |

Both need WSL2 first. Doing both is normal: Docker for running cases and
reproducing someone else's results, a native build for day-to-day work on the
code.

## Setting up WSL2

First check what you already have. On a reasonably current Windows 11 the WSL
runtime is often installed even when no distribution is:

```powershell
wsl --version
wsl --list --verbose
```

If `wsl --version` prints a version and kernel, the runtime is there and you only
need a distribution:

```powershell
wsl --install -d Ubuntu-22.04
```

That downloads roughly 500 MB and opens a shell. **No restart is needed when the
runtime is already installed** - a restart is only required when WSL itself is
being enabled for the first time.

### Two messages that look like problems and are not

**`wsl --status` may say the Windows Subsystem for Linux optional component is
not enabled and WSL1 is unavailable.** That is about WSL**1**, the legacy
translation layer. WSL2 does not use it. If `wsl --version` prints a kernel
version, you are fine.

**`Get-CimInstance Win32_Processor` may report
`VirtualizationFirmwareEnabled: False`.** This is a reporting artefact, not a
BIOS setting. When a hypervisor is already running, Windows itself runs inside
it and cannot read the raw CPU flags. Check this instead:

```powershell
(Get-CimInstance Win32_ComputerSystem).HypervisorPresent
```

`True` means virtualisation is working, whatever the other field says. The
conclusive test is simply that `wsl --list --online` returns a catalogue.

## Two things the first-run setup gets wrong

**The setup prompt can be left unfinished.** `wsl --install` opens a window
asking for a UNIX username. Until you answer it, the distribution is registered
but has no normal user, and every `wsl` command from another terminal blocks
waiting on that prompt. If `wsl` seems to hang, look for the open Ubuntu window.

**Check that a user was actually created.** If setup was skipped or interrupted,
the distribution defaults to root:

```powershell
Get-ItemProperty "HKCU:\SOFTWARE\Microsoft\Windows\CurrentVersion\Lxss\*" |
    Select-Object DistributionName, DefaultUid
```

`DefaultUid` of `0` means root. That causes two problems that surface much later
and look unrelated:

- **OpenMPI refuses to run as root.** Any parallel run fails.
- **Files written into a mounted directory are owned by root**, so the container,
  which runs as uid 1000, cannot write there - and the error you see is a bare
  permission denied several steps from the cause.

Fix it once:

```bash
# inside WSL, as root
adduser yourname
usermod -aG sudo yourname
```

```powershell
# back in PowerShell
ubuntu2204.exe config --default-user yourname
wsl --shutdown
```

## Where to put the repository

This choice has a real cost attached.

| Location | From Windows | Compile speed |
|---|---|---|
| `~/HigFlow` inside WSL | via `\\wsl$\Ubuntu-22.04\home\...` | fast - native ext4 |
| `/mnt/c/dev/HigFlow` | directly, it *is* a Windows folder | slow - every file access crosses a translation layer |

Compiling 176,000 lines across `/mnt/c` is noticeably slower than on the WSL
filesystem, and the gap widens with the number of small files, which is exactly
what a C project is.

**Recommendation:** clone into WSL for anything you compile.

```bash
cd ~
git clone https://github.com/antoniocastelofilho/HigFlow.git
cd HigFlow
```

Windows programs can still reach those files at
`\\wsl$\Ubuntu-22.04\home\yourname\HigFlow` - paste that into Explorer's address
bar, and pin it.

Keeping the repository under `/mnt/c` is a reasonable choice if you mainly edit
documents and rarely compile. Just know why builds feel slow when they do.

## Path A - Docker, for running

Two ways to get Docker, and the choice matters.

**Docker Engine inside WSL** is lighter and free of licence questions. It needs
systemd, which recent Ubuntu images on WSL enable by default:

```bash
cat /etc/wsl.conf
```

should show:

```ini
[boot]
systemd=true
```

If it does not, add it, then `wsl --shutdown` from PowerShell and reopen. Then
follow the Docker Engine instructions in
[the container guide](containers.md#windows). `docker` then works inside WSL,
and from PowerShell as `wsl docker ...`.

**Docker Desktop** integrates with Windows so `docker` works directly in
PowerShell, and gives you a GUI. It is heavier, and requires a paid subscription
for large companies - free for personal use, education and small business.

Either way, from the repository:

```bash
docker build -f containers/Dockerfile -t higflow:latest .
docker run --rm -v "$PWD/cases:/work" higflow:latest case example2d_Newt
```

Measured on sixteen cores: 6 min 30 s to build, of which PETSc is 5 min. Running
the Newtonian case afterwards produced 101 VTK files. The full guide, including
parallel runs and clusters, is in [containers.md](containers.md).

## Path B - native build in WSL, for developing

Inside the Ubuntu shell, follow the Linux instructions. In short:

```bash
sudo apt update
# install the dependencies - see the Linux guide
source varsrc
cd higtree && make DIM=2 && make DIM=3
cd ../higflow && make DIM=2
cd example2d_Newt && make && make run
```

Read [`varsrc`](../../varsrc) before sourcing it: the committed `PETSC_DIR` does
not match where the bundled installer puts PETSc, so it needs correcting to your
actual installation.

One thing that bites specifically here: **the example Makefiles compile with
`gcc`, not `mpicc`**, and take their MPI include path from PETSc's
`PETSC_CC_INCLUDES`. On Ubuntu, OpenMPI's headers are under
`/usr/lib/x86_64-linux-gnu/openmpi/include`, so if PETSc was configured against a
system MPI without that path recorded, every example fails with:

```
pdomain.h:5:10: fatal error: mpi.h: No such file or directory
```

The workaround, until the Makefiles use `mpicc`:

```bash
export CPATH=/usr/lib/x86_64-linux-gnu/openmpi/include
export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib
```

## Viewing results in ParaView on Windows

Install ParaView on **Windows**, not inside WSL. A Windows ParaView is faster,
uses your GPU properly, and opens WSL files without any copying:

1. Install ParaView from [paraview.org](https://www.paraview.org/download/)
2. File → Open
3. Paste `\\wsl$\Ubuntu-22.04\home\yourname\HigFlow\cases\example2d_Newt\VTKS`
4. ParaView groups the numbered `.vtk` files into a single time series - select
   the group, not an individual file, so you can animate

The VTK files HigFlow writes are ASCII `UNSTRUCTURED_GRID`, which both ParaView
and VisIt read directly.

## Editing code from Windows

**VS Code** is the smoothest option. Install the *WSL* extension, then from the
Ubuntu shell:

```bash
code .
```

The editor window runs on Windows; the language server, compiler, terminal and
debugger all run inside Linux. Breakpoints in gdb work.

Any Windows editor can also open `\\wsl$\...` directly, at the cost of the same
filesystem translation that slows compilation.

**Set your line endings before committing anything.** Git for Windows defaults
to `core.autocrlf=true`, which checks shell scripts out with CRLF - and a CRLF
shebang fails on Linux with a message that names the interpreter, not the cause:

```
bash: ./install_higflow_ubuntu22: /bin/bash^M: bad interpreter
```

The repository's `.gitattributes` forces LF on scripts, so a fresh clone is fine.
For a clone made before that existed:

```bash
git config core.autocrlf input
git rm --cached -r . && git reset --hard
```

## Disk and memory

**Disk.** `df -h /` inside WSL reports the virtual disk's maximum size, often
around 1 TB, which is not how much you have. The real limit is free space on
`C:`. Check the honest number:

```bash
df -h /mnt/c
```

Budget roughly 10 GB: the Ubuntu image, the PETSc build, HigFlow, and the
container image at about 1.3 GB. Simulation output is separate and can be large
- the 3D lid-driven case wrote 509 MB of VTK before being stopped.

WSL's virtual disk grows but does not shrink on its own. To reclaim space after
deleting files inside it, run `Optimize-VHD` from an elevated PowerShell, or
`wsl --manage Ubuntu-22.04 --set-sparse true`.

**Memory.** WSL2 takes up to half your RAM by default. To change it, create
`%UserProfile%\.wslconfig`:

```ini
[wsl2]
memory=8GB
processors=8
```

Then `wsl --shutdown` and reopen. Give it at least 4 GB - the PETSc build is
killed by the out-of-memory reaper below that.

## Troubleshooting

**`wsl` hangs with no output.** The first-run setup prompt is open in another
window waiting for a username. Answer it.

**`WslRegisterDistribution failed with error: 0x80370102`.** Virtualisation is
disabled in firmware. Enable Intel VT-x or AMD-V in the BIOS. Distinguish this
from the WMI reporting artefact described earlier: this one is a real failure to
start.

**Docker commands say the daemon is not running.** Inside WSL,
`sudo systemctl start docker`. If systemd is not pid 1 (`ps -p 1 -o comm=`), the
`[boot] systemd=true` setting is missing from `/etc/wsl.conf`.

**Compiling is very slow.** The repository is under `/mnt/c`. Clone into `~`
instead.

**`mpirun` refuses to run as root.** No normal user was created - see
[above](#two-things-the-first-run-setup-gets-wrong).

**Antivirus makes builds crawl.** Real-time scanning inspects every object file.
Exclude the WSL virtual disk directory, under
`%LocalAppData%\wsl`, or `%LocalAppData%\Packages\CanonicalGroupLimited*`.

**A container writes files you cannot delete from Windows.** They are owned by a
uid that does not map to your Windows account. From inside WSL:
`sudo chown -R "$USER:$USER" cases/`.

**Everything worked yesterday and today WSL will not start.** `wsl --shutdown`,
then reopen. A Windows update that replaces the kernel occasionally leaves the
running instance in a bad state.

## MSYS2 and Cygwin

**Not supported, and not recommended to attempt.** MSYS2 gives you a POSIX-ish
toolchain, but the blockers listed at the top do not move: no OpenMPI, no
libfyaml, no libnuma. You would be porting three dependencies before compiling a
line of HigFlow, and maintaining those ports afterwards.

This is stated plainly so nobody spends a weekend discovering it.
