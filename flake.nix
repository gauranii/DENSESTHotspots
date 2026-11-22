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

      # ----- Unit tests -----
      # Shell functions that run the unit tests of the directory they are
      # called in: testthat over tests/testthat/ (unit_r).
      unitTests = pkgs: ''
        export OPENBLAS_NUM_THREADS=1
        export PATH=${rEnv pkgs}/bin:$PATH
        unit_r() { Rscript tests/testthat.R; }
        # Unit tests in other languages are added here, one function each.
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
            # Commands in other languages are listed here.
            ""
            "Other commands:"
            (row "test" "unit tests of the functions in R/")
            (row "versions" "print the software versions")
            (row "help" "print this list (also: nix run . -- --help)")
            "  nix flake check         build every output, compare it with the paper, run the unit tests"
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
            pkgs.python3
          ];
          OPENBLAS_NUM_THREADS = 1;
        };
      });

      apps = forEachSystem (
        pkgs:
        lib.mapAttrs (name: script: app pkgs name [ script ]) scripts
        // {
          all = app pkgs "all" (lib.attrValues scripts);
          versions = {
            type = "app";
            program = "${pkgs.writeShellScript "zinull-versions" ''
              export PATH=${rEnv pkgs}/bin:${pkgs.python3}/bin:$PATH
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
            ''}";
            meta.description = "Run the unit tests of the R code";
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
      # versions differ from README.md; cell by cell results are in the output.
      checks = forEachSystem (pkgs: {
        outputs =
          pkgs.runCommand "zinull-outputs"
            {
              nativeBuildInputs = [
                (rEnv pkgs)
                pkgs.python3
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
      });

      formatter = forEachSystem (pkgs: pkgs.nixfmt-rfc-style);
    };
}
