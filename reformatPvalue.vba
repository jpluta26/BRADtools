'''macro to reformat p-value to publication format
Sub reformatPvalue()

    Dim cell As Range
    Dim inval As String
    Dim parts() As String
    Dim mantissa As String
    Dim exponent As String
    Dim output As String
    Dim expStart As Long

    For Each cell In Selection

        If Not IsEmpty(cell.Value) Then

            inval = Format(cell.Value, "0.00E+00")

            parts = Split(inval, "E")

            mantissa = parts(0)
            exponent = parts(1)

            If Left(exponent, 1) = "+" Then
                exponent = Mid(exponent, 2)
            End If

            If Left(exponent, 1) = "-" Then
                Do While Len(exponent) > 2 And Mid(exponent, 2, 1) = "0"
                    exponent = "-" & Mid(exponent, 3)
                Loop
            Else
                Do While Len(exponent) > 1 And Left(exponent, 1) = "0"
                    exponent = Mid(exponent, 2)
                Loop
            End If

            output = mantissa & " " & ChrW(215) & " 10" & exponent

            cell.Value = output

            expStart = Len(mantissa & " " & ChrW(215) & " 10") + 1

            cell.Characters(Start:=expStart, Length:=Len(exponent)).Font.Superscript = True

        End If

    Next cell

End Sub


