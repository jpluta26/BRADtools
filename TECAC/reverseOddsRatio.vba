
'pluta 12/6/19
'my first attempt at vb

'this is intended for use on the 'working tables' spreadsheet
'reverse odds ratio to align everything to risk allele
'selection a block of stats - include the 4 columns of MAF, blank, OR, OR.ci
'works for an arbitrary number of rows
Sub reverseOR()
 'vb is weird about variable reassignment, probably a better way to do this
    Dim count As Integer
    Dim MAF As Double
    Dim OddsRatio As Double
    Dim newOddsRatio As Double
    Dim StrIn As String
    Dim StrOut() As String
    Dim CIhi As Double
    Dim CIlo As Double
    Dim se1 As Double
    Dim se2 As Double
    Dim se As Double
    Dim outCI(5) As String
  

    Dim rng As Range
    Dim row As Range
    Dim cell As Range
    Set rng = Selection
    

    For Each row In rng.rows
        count = 0
        For Each cell In row.cells
  'the first cell is MAF, invert this
            count = count + 1
            If count = 1 Then
                MAF = 1 - cell.Value
                cell.Value = MAF
   'flip odds ratio
            ElseIf count = 3 Then
                OddsRatio = cell.Value
                newOddsRatio = Exp(-Log(OddsRatio))
                cell.Value = Round(newOddsRatio, 2)
   'derive standard error from confidence interval and apply to inverted OR to get the new confidence interval
   'replace values in existing cells
           ElseIf count = 4 Then
                StrIn = cell.Value
                StrIn = Replace(StrIn, ")", "")
                StrIn = Replace(StrIn, "(", "")
                StrOut = Split(StrIn, ",")
                CIlo = StrOut(0)
                CIhi = StrOut(1)
                se1 = Abs((Log(CIlo) - Log(OddsRatio)) / 1.96)
                se2 = Abs((Log(CIhi) - Log(OddsRatio)) / 1.96)
                se = (se1 + se2) / 2

                outCI(0) = "("
                outCI(1) = Round(Exp(Log(newOddsRatio) - se * 1.96), 2)
                outCI(2) = ", "
                outCI(3) = Round(Exp(Log(newOddsRatio) + se * 1.96), 2)
                outCI(4) = ")"

                cell.Value = Join(outCI, "")
            
            End If
        Next cell
    Next row
End Sub
