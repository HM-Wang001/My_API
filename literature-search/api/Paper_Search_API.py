# api/test.py
import json

def handler(request):
    return {
        'statusCode': 200,
        'headers': {'Content-Type': 'application/json'},
        'body': json.dumps({'message': 'Vercel is working!', 'query': request.get('queryStringParameters', {})})
    }
