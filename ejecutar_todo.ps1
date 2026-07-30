# ============================================================================
#  Ejecucion individual de todos los casos, con guardado incremental.
#  Cada caso corre solo, con todos los hilos, y guarda su propio JSON.
#  Si algo se interrumpe, lo ya calculado no se pierde.
#
#  Uso:   .\ejecutar_todo.ps1
#         .\ejecutar_todo.ps1 -Hilos 16 -Carpeta otros_resultados
# ============================================================================
param(
    [int]$Hilos = 32,
    [string]$Carpeta = "resultados"
)

New-Item -ItemType Directory -Force -Path $Carpeta | Out-Null

# Presupuesto por caso, en segundos. Se asigna segun donde esta la frontera:
#  - casos ya certificados: poco tiempo, solo para reconfirmar
#  - frontera (R=5 z=4, R=4 z=8): mucho tiempo, son los que pueden cerrar
#  - resto: tiempo moderado, solo mejoran cotas superiores
$Presupuesto = @{
    "4-1" =   60; "4-2" =  120; "4-3" =  600; "4-4" = 1800; "4-5" = 14400
    "4-6" = 1800; "4-7" = 1800; "4-8" = 1800; "4-9" = 1800; "4-10" = 1800
    "8-1" =   60; "8-2" =  300; "8-3" = 1800; "8-4" = 14400; "8-5" = 1800
    "8-6" = 1800; "8-7" = 1800; "8-8" = 1800; "8-9" = 1800; "8-10" = 1800
}

$total = 0
foreach ($v in $Presupuesto.Values) { $total += $v }
Write-Host "Presupuesto total: $([math]::Round($total/3600,1)) horas (peor caso)" -ForegroundColor Cyan
Write-Host "Hilos por caso: $Hilos`n" -ForegroundColor Cyan

$inicio = Get-Date
foreach ($z in 4, 8) {
    foreach ($R in 1..10) {
        $clave = "$z-$R"
        $t = $Presupuesto[$clave]
        $salida = Join-Path $Carpeta "z$($z)_R$($R).json"

        if (Test-Path $salida) {
            Write-Host "[omitido] z=$z R=$R  (ya existe $salida)" -ForegroundColor DarkGray
            continue
        }

        Write-Host "`n=== z=$z  R=$R  (limite $t s) ===" -ForegroundColor Yellow
        py src\cpsat_keccak.py --rondas $R --z $z --tiempo $t --hilos $Hilos --salida $salida

        $transcurrido = (Get-Date) - $inicio
        Write-Host "  tiempo acumulado: $([math]::Round($transcurrido.TotalMinutes,1)) min" -ForegroundColor DarkGray
    }
}

Write-Host "`n=== RESUMEN FINAL ===" -ForegroundColor Green
py src\fusionar.py $Carpeta --salida resultados_final.json
