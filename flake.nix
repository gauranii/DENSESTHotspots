{
  description = "Code for Gauran et al. (2018), Biometrics 74:458-471: empirical null estimation with zero-inflated discrete mixtures";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/f3aa7d3d72c2bc66e0f30f4262d7ec2e33c9eab2";

  outputs =
    { self, nixpkgs }:
    let
      lib = nixpkgs.lib;
      systems = [
        "aarch64-darwin"
        "x86_64-darwin"
        "aarch64-linux"
        "x86_64-linux"
      ];
      forEachSystem = f: lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});

      # R with the packages the scripts load (ggplot2) and the package that
      # runs the unit tests (testthat).
      rEnv =
        pkgs:
        pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            ggplot2
            testthat
          ];
        };

      # Python with the packages of python/zinull; numba compiles the random
      # number generator and the loops that follow R's order of operations.
      pyEnv =
        pkgs:
        pkgs.python3.withPackages (
          ps: with ps; [
            numpy
            pandas
            matplotlib
            numba
            pytest
          ]
        );

      # One command per table and figure of the paper (filled in below).
      scripts = import ./scripts.nix;

      # Run the scripts from the flake source; results/ and figures/ are
      # written to the directory the command is started in.
      runner =
        pkgs: name: names:
        pkgs.writeShellScript "zinull-${name}" ''
          set -euo pipefail
          export ZINULL_OUTPUT="''${ZINULL_OUTPUT:-$PWD}"
          export OPENBLAS_NUM_THREADS=1
          # The protein domain counts are not in the repository: use
          # data/realdata.txt of the start directory when it is there.
          if [ -z "''${ZINULL_DATA:-}" ] && [ -f "$PWD/data/realdata.txt" ]; then
            export ZINULL_DATA="$PWD/data/realdata.txt"
          fi
          cd ${self}
          for script in ${lib.concatStringsSep " " names}; do
            ${rEnv pkgs}/bin/Rscript "scripts/$script.R" "$@"
          done
        '';
      app = pkgs: name: names: {
        type = "app";
        program = "${runner pkgs name names}";
        meta.description = "Build ${name} of the paper";
      };

      # The same commands with the Python scripts, as py-<name>.
      pyRunner =
        pkgs: name: names:
        pkgs.writeShellScript "zinull-py-${name}" ''
          set -euo pipefail
          export ZINULL_OUTPUT="''${ZINULL_OUTPUT:-$PWD}"
          export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
          if [ -z "''${ZINULL_DATA:-}" ] && [ -f "$PWD/data/realdata.txt" ]; then
            export ZINULL_DATA="$PWD/data/realdata.txt"
          fi
          export NUMBA_CACHE_DIR="''${NUMBA_CACHE_DIR:-''${TMPDIR:-/tmp}/zinull-numba}"
          export PYTHONDONTWRITEBYTECODE=1
          cd ${self}
          for script in ${lib.concatStringsSep " " names}; do
            ${pyEnv pkgs}/bin/python "python/scripts/$script.py" "$@"
          done
        '';
      pyApp = pkgs: name: names: {
        type = "app";
        program = "${pyRunner pkgs name names}";
        meta.description = "Build ${name} of the paper with the Python scripts";
      };

      # ----- Unit tests -----
      # Shell functions that run the unit tests of the directory they are
      # called in: testthat over tests/testthat/ (unit_r) and pytest over
      # python/tests/ (unit_python).
      unitTests = pkgs: ''
        export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
        export NUMBA_CACHE_DIR="''${NUMBA_CACHE_DIR:-''${TMPDIR:-/tmp}/zinull-numba}"
        export PYTHONDONTWRITEBYTECODE=1
        export PATH=${rEnv pkgs}/bin:${pyEnv pkgs}/bin:$PATH
        unit_r() { Rscript tests/testthat.R; }
        unit_python() { python -m pytest -p no:cacheprovider python/tests; }
      '';

      # ----- Help -----
      # The text `nix run .#help` prints: every command of the flake.
      helpText =
        let
          row =
            command: text:
            "  nix run .#${command}${
                lib.fixedWidthString (lib.max 1 (14 - builtins.stringLength command)) " " ""
              }${text}";
          rows = lib.mapAttrsToList (name: script: row name "scripts/${script}.R") scripts;
        in
        lib.concatStringsSep "\n" (
          [
            "Code of Gauran et al. (2018), Biometrics 74:458-471: empirical null estimation with zero-inflated discrete mixtures"
            ""
            "One command per table and figure of the paper (R):"
          ]
          ++ rows
          ++ [
            (row "all" "every table and figure")
            ""
            "The same commands with the Python scripts (python/scripts/), an independent check of the R numbers:"
            "  nix run .#py-<command>   for example nix run .#py-table1, nix run .#py-all"
            ""
            "Other commands:"
            (row "test" "unit tests of the functions in R/ and python/zinull/")
            (row "versions" "print the software versions")
            (row "help" "print this list (also: nix run . -- --help)")
            "  nix flake check         build every output, compare it with the paper and with Python, run the unit tests"
            "  nix develop             a shell with the same R and Python"
            ""
            "Results are written to results/ and figures/ of the current directory (or of ZINULL_OUTPUT)."
            "Every table and figure command takes the options listed in the header of its script, after --."
            "The simulation tables read the archived replications unless --rerun is given, for example"
            "  nix run .#table1 -- --rerun --runs=1:10 --cores=4"
            "The outputs of the protein domain data need data/realdata.txt in the current directory (or ZINULL_DATA)."
          ]
        );
      help = pkgs: {
        type = "app";
        program = "${pkgs.writeShellScript "zinull-help" ''
          printf '%s\n' "${helpText}"
        ''}";
        meta.description = "List every command";
      };
    in
    {
      devShells = forEachSystem (pkgs: {
        default = pkgs.mkShell {
          packages = [
            (rEnv pkgs)
            (pyEnv pkgs)
          ];
          OPENBLAS_NUM_THREADS = 1;
          VECLIB_MAXIMUM_THREADS = 1;
        };
      });

      apps = forEachSystem (
        pkgs:
        lib.mapAttrs (name: script: app pkgs name [ script ]) scripts
        // lib.mapAttrs' (name: script: lib.nameValuePair "py-${name}" (pyApp pkgs name [ script ])) scripts
        // {
          all = app pkgs "all" (lib.attrValues scripts);
          py-all = pyApp pkgs "all" (lib.attrValues scripts);
          versions = {
            type = "app";
            program = "${pkgs.writeShellScript "zinull-versions" ''
              export PATH=${rEnv pkgs}/bin:${pyEnv pkgs}/bin:$PATH
              cd ${self}
              exec bash checks/versions.sh "$@"
            ''}";
            meta.description = "Print the software versions";
          };
          # The unit tests of the current directory when it is a copy of the
          # repository, of the flake source otherwise.
          test = {
            type = "app";
            program = "${pkgs.writeShellScript "zinull-test" ''
              set -euo pipefail
              ${unitTests pkgs}
              [ -f tests/testthat.R ] || cd ${self}
              unit_r
              unit_python
            ''}";
            meta.description = "Run the unit tests of the R and the Python code";
          };
          # `nix run .#help`, `nix run . -- --help` and `nix run .` print the
          # list of commands.
          help = help pkgs;
          default = help pkgs;
        }
      );

      # Build every table and figure with the default settings and compare
      # each with the paper; the result is the comparison (checks/RESULTS.md).
      # The check fails when an output cannot be built or the software
      # versions differ from docs/DETAILS.md; cell by cell results are in the output.
      checks = forEachSystem (pkgs: {
        outputs =
          pkgs.runCommand "zinull-outputs"
            {
              nativeBuildInputs = [
                (rEnv pkgs)
                (pyEnv pkgs)
              ];
            }
            ''
              export HOME=$TMPDIR ZINULL_OUTPUT=$TMPDIR
              cd ${self}
              bash checks/run_checks.sh $TMPDIR/RESULTS.md || test -s $TMPDIR/RESULTS.md
              bash checks/versions.sh --check
              cp $TMPDIR/RESULTS.md $out
            '';

        # The unit tests of every function in R/ (testthat, tests/testthat/).
        unit-r = pkgs.runCommand "zinull-unit-r" { } ''
          export HOME=$TMPDIR
          ${unitTests pkgs}
          cd ${self}
          unit_r | tee $TMPDIR/unit-r.txt
          cp $TMPDIR/unit-r.txt $out
        '';

        # The unit tests of every function in python/zinull/ (pytest,
        # python/tests/); the values of R they compare with are stored in
        # python/tests/r_values.json.
        unit-python = pkgs.runCommand "zinull-unit-python" { nativeBuildInputs = [ (pyEnv pkgs) ]; } ''
          export HOME=$TMPDIR NUMBA_CACHE_DIR=$TMPDIR PYTHONDONTWRITEBYTECODE=1
          export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
          set -o pipefail
          cd ${self}
          python -m pytest -p no:cacheprovider python/tests | tee $TMPDIR/unit-python.txt
          cp $TMPDIR/unit-python.txt $out
        '';

        # Build every table and figure with the R and with the Python
        # scripts and compare the two; the result is the comparison
        # (checks/PYTHON.md, with the stored comparison of the reruns).
        python =
          pkgs.runCommand "zinull-python"
            {
              nativeBuildInputs = [
                (rEnv pkgs)
                (pyEnv pkgs)
              ];
            }
            ''
              export HOME=$TMPDIR NUMBA_CACHE_DIR=$TMPDIR PYTHONDONTWRITEBYTECODE=1
              export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 CORES=4
              cd ${self}
              bash checks/check_python.sh archive $TMPDIR/work $TMPDIR/PYTHON.md $TMPDIR \
                || { cat $TMPDIR/PYTHON.md; exit 1; }
              cp $TMPDIR/PYTHON.md $out
            '';
      });

      formatter = forEachSystem (pkgs: pkgs.nixfmt-rfc-style);
    };
}
