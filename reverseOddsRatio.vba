
'given a MAF, odds ratio, and confidence interval, invert all values
'this is the same as flipping the reference allele
'use this code to align all variants to the RISK allele
Sub reverseOR()
'VB is weird about variable assignment... probably a better way to do this
 
    Dim cell As Object
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
  
    count = 0

    'select the four column range containing the values
    For Each cell In Selection
    'first cell if RAF
        count = count + 1
        If count = 1 Then
           MAF = 1 - cell.Value
           cell.Value = MAF
    'second cell is OR
        ElseIf count = 3 Then
            OddsRatio = cell.Value
            newOddsRatio = Exp(-Log(OddsRatio))
            cell.Value = Round(newOddsRatio, 2)
    'fourth cell is confidence interval of OR
        ElseIf count = 4 Then
    'derive standard error of logOR from the confidence interval
            StrIn = cell.Value
            StrIn = Replace(StrIn, ")", "")
            StrIn = Replace(StrIn, "(", "")
            StrOut = Split(StrIn, ",")
            CIlo = StrOut(0)
            CIhi = StrOut(1)
            se1 = Abs((Log(CIlo) - Log(OddsRatio)) / 1.96)
            se2 = Abs((Log(CIhi) - Log(OddsRatio)) / 1.96)
            se = (se1 + se2) / 2

    'then apply to the inverted OR for the inverted CI
            outCI(0) = "("
            outCI(1) = Round(Exp(Log(newOddsRatio) - se * 1.96), 2)
            outCI(2) = ", "
            outCI(3) = Round(Exp(Log(newOddsRatio) + se * 1.96), 2)
            outCI(4) = ")"

    'write new value in place
            cell.Value = Join(outCI, "")
            
        End If
     Next cell

End Sub


