module() {
  case "$1" in
    -t) printf '%s\n' "/opt/apps/eb/modules/all:" "Python/2.7.18-GCCcore-12.2.0" "Python/3.10.8-GCCcore-12.2.0" "Python/3.10.8-GCCcore-12.2.0-bare" "Python/3.11.3-GCCcore-12.3.0 (D)" "Python/3.12.3-GCCcore-13.3.0(default)" >&2 ;;
    load) echo "fake module load [$2]" >&2; case "$2" in Python/3.11.3-GCCcore-12.3.0|Python/3.10.8-GCCcore-12.2.0) export PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/newpy:$PATH";; *) return 1;; esac ;;
    unload) ;;
  esac
}
